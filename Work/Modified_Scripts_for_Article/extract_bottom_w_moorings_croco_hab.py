'''
Model velocity at the 7 RREX moorings, at the SAME HEIGHT ABOVE BOTTOM (hab) as the deepest current
meter of each mooring.  Uses only Python_Modules_p3 (load, var, R_tools, R_tools_fort).

Which hab?
    observations : hab_obs = DEPTH_REF - z_sensor
                   z_sensor  = median depth of the deepest bin of the mooring (RREX_<name>_CURR_hourly.nc)
                   DEPTH_REF = surveyed water depth at the mooring
    model        : z = -h_model + hab_obs, i.e. the same height above the MODEL sea floor at the grid
                   point of the mooring (the model floor can differ from the real one by tens of
                   metres: equal hab is not equal depth).

How
    The observed bottom slope is computed from the SRTM15 bathymetry (TOPO15_FILE) at several smoothing
    lengths (SLOPE_LENGTHS_KM); OBS_SLOPE_KM is the length used for the observed w_u and for the model
    slope smoothing.  Mooring positions come from the mooring files themselves.

    For each mooring a small window (+-HALF grid points) around its nearest model point is read, so that
    only a few columns are loaded per snapshot.  u, v (rho points), the true vertical velocity w
    (var('w'): continuity, toolsF.get_wvlcty) and the depths (tools.get_depths) of the window are
    interpolated linearly onto HAB_GRID.  The bottom slope is taken from the same window, at grid scale
    and smoothed over SMOOTH_KM (same length as the observed slope, OBS_SLOPE_KM).

Output (netCDF, dims snapshot, mooring, hab)
    u, v, w            (snapshot, mooring, hab)  [m s-1]
    wu_obs_hab, wu_s_obs_hab, w_obs_hab  (snapshot, mooring)  [m day-1], at hab_obs:
          wu    = -u.grad(h)             grid-scale slope
          wu_s  = -u.grad(h_smoothed)    slope smoothed over SMOOTH_KM, like the observations
          w     = model true vertical velocity
    dhdx, dhdy  (mooring)            model slope at grid scale; dhdx_s, dhdy_s: smoothed over SMOOTH_KM
    dhdx_L, dhdy_L (length, mooring) model slope smoothed over each SLOPE_LENGTHS_KM (0 = grid scale), like the obs
    wu_L_obs_hab (snapshot, mooring, length) [m day-1] = -u.grad(h smoothed over L) at hab_obs
    h, hab_obs, z_sensor, depth_ref, i_model, j_model  (mooring)
    slope_length_km (length); D_topo15 (mooring); dDdx_topo15, dDdy_topo15 (length, mooring)   
          = bathymetry slopes from SRTM15 smoothed over each length (east, north components)
'''
import os, sys
import numpy as np
from netCDF4 import Dataset
from scipy.ndimage import uniform_filter

# folder of Python_Modules_p3 if it is not already on PYTHONPATH (e.g. on datarmor), else None
MODULES_DIR = os.environ.get('PYTHON_MODULES_P3', None)
if MODULES_DIR: sys.path.append(MODULES_DIR)

from Modules import *            # load
from Modules_gula import *       # var
import R_tools as tools
import R_tools_fort as toolsF    # (imported by Modules too; kept explicit)

# ------------------------------------------------------------------ configuration
# simulation: first command-line argument, else the 12-hourly run.  Two runs are available:
#   rrexnumsb200_RSUP5_NOFILT_tracer                  (output every 12 h)
#   rrexnumsb200_RSUP5_NOFILT_tracer_article_hourly   (output every 1 h)
SIMUL     = sys.argv[1] if len(sys.argv) > 1 else 'rrexnumsb200_RSUP5_NOFILT_tracer'
# period: from the first snapshot below (spin-up skipped) to the last history file found in the output folder of the run.
# The time step between 2 snapshots is simul.dt [s] (read from the files).
FIRST_SNAPSHOT = {'rrexnumsb200_RSUP5_NOFILT_tracer': 22,                  # rrexnum200_his.00022.nc (avoids the spin-up)
                  'rrexnumsb200_RSUP5_NOFILT_tracer_article_hourly': 0}    # rrexnum200_his_1h.00000.nc: no spin-up, continues the
                                                                           # 12 h run (2nd time step of rrexnum200_his.00060.nc + 1 h)
FIRST_DEFAULT  = 0                                           # first snapshot for a run that is not in FIRST_SNAPSHOT
SNAPSHOT_RANGE = None                                        # (first, last) to override the period
DT_OVERRIDE = None                                           # output interval [s] if simul.dt is not available (e.g. 3600.)

MOORINGS  = ['IRW', 'IRM', 'IRE', 'RRT', 'ICW', 'ICM', 'ICE']
# in-situ files: RREX_<name>_CURR_hourly.nc, topo15_NorthAtl.nc
#   datarmor : '/home/datawork-lops-rrex/nschifan/Data_in_situ_Rene/'
DATA_DIR  = './data/'
FILE_OUT  = None                                             # None: ./rrex_w_moorings_hab.nc (12 h run) or ./rrex_w_moorings_hab_<dt>h.nc (other output interval)
# surveyed water depth at the mooring [m]
DEPTH_REF = np.array([2160., 1570., 1580., 1528., 2162., 2151., 2379.])

TOPO15_FILE = DATA_DIR + 'topo15_NorthAtl.nc'                # SRTM15+ (lat, lon, z<0)
SLOPE_LENGTHS_KM = [0., 2., 5., 10., 20.]                    # smoothing lengths computed for the observed slope (0 = none)
OBS_SLOPE_KM = 10.                                           # length used for w_u (must be in SLOPE_LENGTHS_KM)
EARTH_R   = 6371e3

HALF      = 20                                               # window half-width [grid points]; must exceed half the largest smoothing box (checked)
SMOOTH_KM = OBS_SLOPE_KM                                     # length of the model slope saved as dhdx_s / wu_s (all SLOPE_LENGTHS_KM are saved too)
HAB_GRID  = np.arange(0., 302., 2.)                          # heights above bottom saved [m]
NZ        = 200                                              # sigma levels in the files
SEC_PER_DAY = 86400.


# ------------------------------------------------------------------ observations
def deepest_sensor_depth(name):
    """Median depth [m] of the deepest bin (largest median depth) of a mooring."""
    with Dataset(DATA_DIR + 'RREX_%s_CURR_hourly.nc' % name) as ds:
        d = np.ma.filled(ds['DEPTH'][:].astype(float), np.nan)
    d[np.abs(d) > 1e5] = np.nan
    return float(np.nanmax(np.nanmedian(d, axis=0)))


def topo15_slopes(file, lon_m, lat_m, lengths_km, half_deg=0.6):
    """Bathymetry slopes at the moorings from a regular lat/lon topography (SRTM15, ~0.46 km).

    For each smoothing length L: D = -z is averaged with a square boxcar of L x L km (sizes in grid points
    from the local metric: dy = R dlat, dx = R cos(lat) dlon), then differentiated with centred differences,
    dD/dx = (1/(R cos lat)) dD/dlon (eastward), dD/dy = (1/R) dD/dlat (northward), and interpolated bilinearly
    to the mooring.  L = 0 : no smoothing.
    Returns D (nm,), dDdx (nL, nm), dDdy (nL, nm)."""
    from scipy.interpolate import RegularGridInterpolator
    nm, nL = len(lon_m), len(lengths_km)
    D, dDdx, dDdy = np.full(nm, np.nan), np.full((nL, nm), np.nan), np.full((nL, nm), np.nan)
    with Dataset(file) as nc:
        lat, lon = np.asarray(nc['lat'][:], float), np.asarray(nc['lon'][:], float)
        for m, (lo, la) in enumerate(zip(lon_m, lat_m)):
            j = np.flatnonzero((lat >= la - half_deg) & (lat <= la + half_deg))
            i = np.flatnonzero((lon >= lo - half_deg) & (lon <= lo + half_deg))
            z = np.ma.filled(nc['z'][j[0]:j[-1] + 1, i[0]:i[-1] + 1].astype(float), np.nan)
            la_w, lo_w = lat[j[0]:j[-1] + 1], lon[i[0]:i[-1] + 1]
            Dw = -z                                                  # positive depth
            dlat, dlon = np.deg2rad(np.mean(np.diff(la_w))), np.deg2rad(np.mean(np.diff(lo_w)))
            dy = EARTH_R * dlat                                      # [m]
            dx = EARTH_R * np.cos(np.deg2rad(la_w)) * dlon           # [m], one value per row
            D[m] = RegularGridInterpolator((la_w, lo_w), Dw)([la, lo])[0]
            for k, L in enumerate(lengths_km):
                if L > 0:
                    n = (max(1, int(round(L * 1e3 / dy))), max(1, int(round(L * 1e3 / np.interp(la, la_w, dx)))))
                    Ds = uniform_filter(Dw, size=n, mode='nearest')
                else:
                    Ds = Dw
                gy = np.gradient(Ds, axis=0) / dy
                gx = np.gradient(Ds, axis=1) / dx[:, None]
                dDdx[k, m] = RegularGridInterpolator((la_w, lo_w), gx)([la, lo])[0]
                dDdy[k, m] = RegularGridInterpolator((la_w, lo_w), gy)([la, lo])[0]
    return D, dDdx, dDdy


def mooring_position(name):
    """(lon, lat) of a mooring from its data file."""
    with Dataset(DATA_DIR + 'RREX_%s_CURR_hourly.nc' % name) as ds:
        return float(np.ravel(ds['LONGITUDE'][:])[0]), float(np.ravel(ds['LATITUDE'][:])[0])


# ------------------------------------------------------------------ model windows
def window_around(x, y, lon, lat, half):
    """Nearest model point to (lon, lat) and the window of +-half points around it.
    Returns (i, j, coordinates string for load(), i_min, j_min)."""
    d2 = (np.cos(np.deg2rad(lat)) * (x - lon)) ** 2 + (y - lat) ** 2
    i, j = np.unravel_index(np.argmin(d2), d2.shape)
    i0, j0 = max(i - half, 0), max(j - half, 0)
    coord = ' [%d,%d,%d,%d,[1,%d,1]] ' % (j0, j + half + 1, i0, i + half + 1, NZ)    # [j range, i range, levels]
    return int(i), int(j), coord, i0, j0

def bathymetry_slopes(h, pm, pn, length_km):
    """dh/dx, dh/dy [-] at grid scale and after a boxcar average of h over length_km."""
    n = (max(1, int(round(length_km * 1e3 * np.nanmean(pm)))), max(1, int(round(length_km * 1e3 * np.nanmean(pn)))))
    h_s = uniform_filter(np.asarray(h, float), size=n, mode='nearest')
    grad = lambda f: (np.gradient(f, axis=0) * pm, np.gradient(f, axis=1) * pn)
    return grad(h), grad(h_s)

def to_hab(f_col, hab_col, hab_grid=HAB_GRID):
    """Linear interpolation of a column (hab increasing) onto hab_grid; NaN outside, no extrapolation."""
    return np.interp(hab_grid, hab_col, f_col, left=np.nan, right=np.nan)

def column_profiles(sim, i, j):
    """u, v, w [m s-1] of the column (i, j) of the window `sim`, on HAB_GRID."""
    z_r, z_w = tools.get_depths(sim)
    u = tools.u2rho(var('u', sim).data)
    v = tools.v2rho(var('v', sim).data)
    w = var('w', sim).data                                   # true w from continuity
    floor = z_w[i, j, 0]
    h_r = z_r[i, j, :] - floor
    h_w = z_w[i, j, :] - floor if w.shape[2] == z_w.shape[2] else h_r      # w on w-levels or rho-levels
    return to_hab(u[i, j, :], h_r), to_hab(v[i, j, :], h_r), to_hab(w[i, j, :], h_w)


def available_snapshots(sim, tfile):
    """Snapshot indices covered by the history files present in the folder of sim.ncfile
    (file number F holds snapshots F .. F+tfile-1).  Missing files inside the period are skipped later."""
    import re
    folder, base = os.path.split(sim.ncfile)
    m = re.match(r'^(.*\.)(\d+)(\.nc)$', base)
    assert m, 'cannot parse the history file name %s' % base
    pat = re.compile('^' + re.escape(m.group(1)) + r'(\d+)' + re.escape(m.group(3)) + '$')
    nums = [int(g.group(1)) for g in (pat.match(f) for f in os.listdir(folder or '.')) if g is not None]
    assert nums, 'no history file %s*%s in %s' % (m.group(1), m.group(3), folder)
    return min(nums), max(nums) + tfile - 1


_NREC = {}
def has_snapshot(sim):
    """True if sim.ncfile exists and holds the record sim.infiletime (the last file of a run can be incomplete)."""
    f = sim.ncfile
    if f not in _NREC:
        if not os.path.exists(f):
            _NREC[f] = 0
        else:
            with Dataset(f) as nc:
                names = [n for n in ('time', 'time_counter', 'ocean_time', 'scrum_time') if n in nc.variables]
                _NREC[f] = len(nc.variables[names[0]]) if names else \
                           max([len(d) for k, d in nc.dimensions.items() if k.startswith('time')] or [10 ** 9])
    return sim.infiletime < _NREC[f]


# ------------------------------------------------------------------ main
def main():
    nm, nh = len(MOORINGS), len(HAB_GRID)
    lon_m, lat_m = np.array([mooring_position(n) for n in MOORINGS]).T
    z_sensor = np.array([deepest_sensor_depth(n) for n in MOORINGS])
    hab_obs = DEPTH_REF - z_sensor

    # --- observed bottom slope from SRTM15
    D15, dDdx15, dDdy15 = topo15_slopes(TOPO15_FILE, lon_m, lat_m, SLOPE_LENGTHS_KM)
    topo = (D15, dDdx15, dDdy15)
    k_obs = SLOPE_LENGTHS_KM.index(OBS_SLOPE_KM)
    print('depth at the mooring [m]: SRTM15 | DEPTH_REF')
    for m, n in enumerate(MOORINGS):
        print('  %-4s lon %8.3f lat %7.3f %8.0f %8.0f' % (n, lon_m[m], lat_m[m], D15[m], DEPTH_REF[m]))
    for name, arr in (('dDdx (east)', dDdx15), ('dDdy (north)', dDdy15)):
        print('\n%s from SRTM15 smoothed over L (used: %g km)\n      ' % (name, OBS_SLOPE_KM) + ''.join('%9s' % ('%gkm' % L) for L in SLOPE_LENGTHS_KM))
        for m, n in enumerate(MOORINGS):
            print('  %-4s' % n + ''.join('%9.4f' % arr[k, m] for k in range(len(SLOPE_LENGTHS_KM))))
    print()

    # --- one small window per mooring
    t_start = SNAPSHOT_RANGE[0] if SNAPSHOT_RANGE else FIRST_SNAPSHOT.get(SIMUL, FIRST_DEFAULT)
    grid = load(SIMUL, time=t_start, output=False)                    # full grid (x, y) to locate the moorings
    dt = float(DT_OVERRIDE or grid.dt)                                # output interval of the run [s] (simul.dt: difference of the last two times of the file)
    assert dt > 0, 'simul.dt = %g: file with a single record, set DT_OVERRIDE (seconds)' % dt
    t_first, t_last = t_start, (SNAPSHOT_RANGE[1] if SNAPSHOT_RANGE else available_snapshots(grid, int(grid.ncname.tfile))[1])
    SNAPSHOTS = list(range(t_first, t_last + 1))
    file_out = FILE_OUT or ('./rrex_w_moorings_hab.nc' if abs(dt - 43200.) < 1 else './rrex_w_moorings_hab_%gh.nc' % (dt / 3600.))
    print('simulation %s: output every %g h, %d snapshots/file -> snapshots %d..%d (%d, %.1f days), written to %s\n'
          % (SIMUL, dt / 3600., int(grid.ncname.tfile), SNAPSHOTS[0], SNAPSHOTS[-1], len(SNAPSHOTS), len(SNAPSHOTS) * dt / 86400., file_out))
    wins, ij, ij_loc = [], [], []
    for lo, la in zip(lon_m, lat_m):
        i, j, coord, i0, j0 = window_around(grid.x, grid.y, lo, la, HALF)
        wins.append(load(simul=grid.simul + coord + format(SNAPSHOTS[0]), output=False))
        ij.append((i, j)); ij_loc.append((i - i0, j - j0))
    ii, jj = np.array(ij).T

    # --- bathymetry: h, slopes
    h_m = np.array([w_.topo[a, b] for w_, (a, b) in zip(wins, ij_loc)])
    # model slope for every smoothing length of the observations (0 km = grid scale)
    for w_ in wins:
        need = int(round(max(SLOPE_LENGTHS_KM) * 1e3 * np.nanmean(w_.pm) / 2.)) + 2
        if need > HALF: print('WARNING: HALF=%d too small for the %g km box (needs >= %d): edge effects' % (HALF, max(SLOPE_LENGTHS_KM), need))
    dhdx_L = np.zeros((len(SLOPE_LENGTHS_KM), nm)); dhdy_L = np.zeros_like(dhdx_L)
    for k, L in enumerate(SLOPE_LENGTHS_KM):
        for m, (w_, (a, b)) in enumerate(zip(wins, ij_loc)):
            gx, gy = bathymetry_slopes(w_.topo, w_.pm, w_.pn, L)[1]    # boxcar over L km, then gradient (L=0: grid scale)
            dhdx_L[k, m], dhdy_L[k, m] = gx[a, b], gy[a, b]
    k_s = SLOPE_LENGTHS_KM.index(SMOOTH_KM)
    k_0 = SLOPE_LENGTHS_KM.index(0.) if 0. in SLOPE_LENGTHS_KM else None
    dhdx_s, dhdy_s = dhdx_L[k_s], dhdy_L[k_s]
    dhdx = dhdx_L[k_0] if k_0 is not None else np.array([bathymetry_slopes(w_.topo, w_.pm, w_.pn, 0)[0][0][a, b] for w_, (a, b) in zip(wins, ij_loc)])
    dhdy = dhdy_L[k_0] if k_0 is not None else np.array([bathymetry_slopes(w_.topo, w_.pm, w_.pn, 0)[0][1][a, b] for w_, (a, b) in zip(wins, ij_loc)])
    for name, arr in (('dhdx (east)', dhdx_L), ('dhdy (north)', dhdy_L)):
        print('model %s, boxcar over L (obs: dDd* above)\n      ' % name + ''.join('%9s' % ('%gkm' % L) for L in SLOPE_LENGTHS_KM))
        for m, n in enumerate(MOORINGS):
            print('  %-4s' % n + ''.join('%9.4f' % arr[k, m] for k in range(len(SLOPE_LENGTHS_KM))))
    print()

    print('%-4s %9s %10s %9s | %8s | %10s' % ('', 'z_sensor', 'depth_ref', 'h_model', 'hab_obs', 'model-real'))
    for k, n in enumerate(MOORINGS):
        print('%-4s %9.1f %10.0f %9.1f | %8.1f | %10.1f' % (n, z_sensor[k], DEPTH_REF[k], h_m[k], hab_obs[k], h_m[k] - DEPTH_REF[k]))
    print('model sampled at hab_obs above the model floor (z = -h_model + hab_obs)\n')

    # --- time loop
    nt = len(SNAPSHOTS)
    U, V, W = (np.full((nt, nm, nh), np.nan, np.float32) for _ in range(3))
    done = np.zeros(nt, bool)
    for n, t in enumerate(SNAPSHOTS):
        for m, (sim, (a, b)) in enumerate(zip(wins, ij_loc)):
            sim.update(time=t)
            if not has_snapshot(sim):
                print('missing or incomplete file for snapshot %d: %s' % (t, sim.ncfile)); break
            U[n, m], V[n, m], W[n, m] = column_profiles(sim, a, b)
        else:
            done[n] = True
        if n % 50 == 0: print('snapshot %d (%d/%d)' % (t, n + 1, nt))

    # --- series at hab_obs
    at_hab = lambda F: np.array([[np.interp(hab_obs[m], HAB_GRID, F[n, m], left=np.nan, right=np.nan)
                                  for m in range(nm)] for n in range(nt)])
    u_h, v_h, w_h = at_hab(U), at_hab(V), at_hab(W)
    w_obs_hab    = w_h * SEC_PER_DAY
    wu_obs_hab   = -(u_h * dhdx   + v_h * dhdy)   * SEC_PER_DAY
    wu_s_obs_hab = -(u_h * dhdx_s + v_h * dhdy_s) * SEC_PER_DAY
    wu_L_obs_hab = -(u_h[:, :, None] * dhdx_L.T[None] + v_h[:, :, None] * dhdy_L.T[None]) * SEC_PER_DAY   # (snapshot, mooring, length)

    # --- save
    with Dataset(file_out, 'w') as nc:
        nc.createDimension('snapshot', nt); nc.createDimension('mooring', nm); nc.createDimension('hab', nh)
        nc.simulation, nc.moorings, nc.smooth_km = SIMUL, ','.join(MOORINGS), SMOOTH_KM
        nc.dt_s = dt                                             # output interval of the run [s]
        nc.comment = 'model sampled at hab_obs = depth_ref - z_sensor above the model sea floor'
        def put(name, dims, arr, units=''):
            v_ = nc.createVariable(name, 'f4', dims); v_[:] = arr; v_.units = units
        put('snapshot', ('snapshot',), np.array(SNAPSHOTS))
        put('time_days', ('snapshot',), np.array(SNAPSHOTS) * dt / 86400., 'days')   # time since the start of the run
        put('hab', ('hab',), HAB_GRID, 'm')
        for name, arr in (('u', U), ('v', V), ('w', W)):
            put(name, ('snapshot', 'mooring', 'hab'), arr, 'm s-1')
        for name, arr in (('dhdx', dhdx), ('dhdy', dhdy), ('dhdx_s', dhdx_s), ('dhdy_s', dhdy_s), ('h', h_m),
                          ('hab_obs', hab_obs), ('z_sensor', z_sensor), ('depth_ref', DEPTH_REF), ('i_model', ii), ('j_model', jj)):
            put(name, ('mooring',), arr)
        nc.obs_slope_km = OBS_SLOPE_KM
        nc.createDimension('length', len(SLOPE_LENGTHS_KM))
        put('slope_length_km', ('length',), SLOPE_LENGTHS_KM, 'km')
        put('lon', ('mooring',), lon_m); put('lat', ('mooring',), lat_m)
        put('D_topo15', ('mooring',), topo[0], 'm')
        put('dDdx_topo15', ('length', 'mooring'), topo[1]); put('dDdy_topo15', ('length', 'mooring'), topo[2])
        put('dhdx_L', ('length', 'mooring'), dhdx_L); put('dhdy_L', ('length', 'mooring'), dhdy_L)   # model slopes, same lengths as obs
        put('wu_L_obs_hab', ('snapshot', 'mooring', 'length'), wu_L_obs_hab, 'm day-1')
        # observed w_u at the sensor from these slopes is computed in the plotting notebook
        for name, arr in (('w_obs_hab', w_obs_hab), ('wu_obs_hab', wu_obs_hab), ('wu_s_obs_hab', wu_s_obs_hab)):
            put(name, ('snapshot', 'mooring'), arr, 'm day-1')
    print('saved', file_out, '| complete snapshots: %d/%d' % (done.sum(), nt))


if __name__ == '__main__':
    main()
