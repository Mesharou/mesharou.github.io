'''
Maps of the time standard deviation (and mean) of the near-bottom vertical velocity of a CROCO run, for Figure 5.
Uses only Python_Modules_p3 (load, var, R_tools).

At the lowest rho-level (and at the heights HAB_EXTRA above the bottom) of every snapshot of the run (all the history files found, spin-up skipped):
    w   [m day-1]  true vertical velocity (var('w'): continuity)
    wu  [m day-1]  -u.grad(h), the velocity that follows the bottom slope (u, v of the lowest level, grid-scale slope of h)
The statistics are accumulated snapshot by snapshot (Welford), so the time series are never stored.

Usage     python extract_bottom_w_std_map.py [simulation]       (default: the 12-hourly run)
Output    ./rrex_bottom_w_std_map[_<dt>h].nc   dimensions xi, eta
          std_w, mean_w, std_wu, mean_wu (lowest rho-level), the same with the suffix _hab<h> for each height h of HAB_EXTRA
          (e.g. std_w_hab50), count (number of snapshots used), h, pm, pn, mask, lon, lat
          attributes: simulation, dt_s, stride, first_snapshot, last_snapshot
Cost      one full-domain 3D velocity computation per snapshot: use SNAPSHOT_STRIDE for the hourly run.
'''
import os, re, sys
import numpy as np
from netCDF4 import Dataset

MODULES_DIR = os.environ.get('PYTHON_MODULES_P3', None)     # folder of Python_Modules_p3 if not already on PYTHONPATH
if MODULES_DIR: sys.path.append(MODULES_DIR)

from Modules import *            # load
from Modules_gula import *       # var
import R_tools as tools

# ------------------------------------------------------------------ configuration
SIMUL = sys.argv[1] if len(sys.argv) > 1 else 'rrexnumsb200_RSUP5_NOFILT_tracer'
#   rrexnumsb200_RSUP5_NOFILT_tracer                  (output every 12 h)
#   rrexnumsb200_RSUP5_NOFILT_tracer_article_hourly   (output every 1 h)
FIRST_SNAPSHOT = {'rrexnumsb200_RSUP5_NOFILT_tracer': 22,                 # rrexnum200_his.00022.nc (avoids the spin-up)
                  'rrexnumsb200_RSUP5_NOFILT_tracer_article_hourly': 0}   # rrexnum200_his_1h.00000.nc
FIRST_DEFAULT = 0
# use one snapshot out of SNAPSHOT_STRIDE (the hourly run has 12 times more snapshots: 12 gives the sampling of the 12-hourly run)
SNAPSHOT_STRIDE = {'rrexnumsb200_RSUP5_NOFILT_tracer': 1,
                   'rrexnumsb200_RSUP5_NOFILT_tracer_article_hourly': 12}
STRIDE_DEFAULT = 1
SNAPSHOT_RANGE = None            # (first, last) to override the period
DT_OVERRIDE = None               # output interval [s] if simul.dt is not available
FILE_OUT = None                  # None: ./rrex_bottom_w_std_map.nc (12 h run) or ./rrex_bottom_w_std_map_<dt>h.nc
HAB_EXTRA = [50., 100.]          # also compute the statistics at these heights above the bottom [m] (linear interpolation in the
                                 # vertical, grid-scale slope of h for wu); [] = lowest model level only
SEC_PER_DAY = 86400.


def available_snapshots(sim, tfile):
    """Snapshot indices covered by the history files present in the folder of sim.ncfile (file F holds F .. F+tfile-1)."""
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


def at_height(F, h_lev, hab):
    """F (nx, ny, nk) linearly interpolated at the height hab above the bottom; h_lev: height of the levels above the bottom
    (nx, ny, nk, increasing).  NaN where hab is outside the water column."""
    k = np.clip((h_lev <= hab).sum(axis=2) - 1, 0, F.shape[2] - 2)[:, :, None]
    h0, h1 = (np.take_along_axis(h_lev, k + dk, 2)[:, :, 0] for dk in (0, 1))
    f0, f1 = (np.take_along_axis(F, k + dk, 2)[:, :, 0] for dk in (0, 1))
    with np.errstate(divide='ignore', invalid='ignore'):
        out = f0 + (hab - h0) / (h1 - h0) * (f1 - f0)
    out[(hab < h_lev[:, :, 0]) | (hab > h_lev[:, :, -1])] = np.nan
    return out


def bottom_fields(sim, habs=()):
    """Dictionary {key: (w, -u.grad(h))} in m day-1, arrays (nx, ny): key None = lowest rho-level, key = hab [m] for each
    height above the bottom of `habs` (linear interpolation of u, v, w on the vertical)."""
    h = np.asarray(sim.topo, float)
    dhdx = np.gradient(h, axis=0) * sim.pm
    dhdy = np.gradient(h, axis=1) * sim.pn
    u = tools.u2rho(var('u', sim).data)
    v = tools.v2rho(var('v', sim).data)
    w = var('w', sim).data                                    # rho-levels: index 0 = lowest level
    out = {None: (w[:, :, 0] * SEC_PER_DAY, -(u[:, :, 0] * dhdx + v[:, :, 0] * dhdy) * SEC_PER_DAY)}
    if len(habs):
        z_r, z_w = tools.get_depths(sim)
        h_r = z_r - z_w[:, :, :1]                             # height of the rho-levels above the bottom
        h_w = h_r if w.shape[2] == z_r.shape[2] else z_w - z_w[:, :, :1]
        for hab in habs:
            uh, vh, wh = at_height(u, h_r, hab), at_height(v, h_r, hab), at_height(w, h_w, hab)
            out[hab] = (wh * SEC_PER_DAY, -(uh * dhdx + vh * dhdy) * SEC_PER_DAY)
    return out


class Welford:
    """Running count / mean / variance of a field, ignoring NaN."""
    def __init__(self, shape):
        self.n = np.zeros(shape); self.mean = np.zeros(shape); self.m2 = np.zeros(shape)
    def add(self, x):
        ok = np.isfinite(x)
        self.n += ok
        delta = np.where(ok, x - self.mean, 0.)
        self.mean += np.divide(delta, self.n, out=np.zeros_like(delta), where=self.n > 0)
        self.m2 += delta * np.where(ok, x - self.mean, 0.)
    def std(self):
        with np.errstate(invalid='ignore', divide='ignore'):
            return np.where(self.n > 1, np.sqrt(self.m2 / self.n), np.nan)       # population std, like np.std
    def mean_(self):
        return np.where(self.n > 0, self.mean, np.nan)


def main():
    t_start = SNAPSHOT_RANGE[0] if SNAPSHOT_RANGE else FIRST_SNAPSHOT.get(SIMUL, FIRST_DEFAULT)
    sim = load(SIMUL, time=t_start, output=False)               # full domain
    dt = float(DT_OVERRIDE or sim.dt)
    assert dt > 0, 'simul.dt = %g: set DT_OVERRIDE (seconds)' % dt
    t_last = SNAPSHOT_RANGE[1] if SNAPSHOT_RANGE else available_snapshots(sim, int(sim.ncname.tfile))[1]
    stride = SNAPSHOT_STRIDE.get(SIMUL, STRIDE_DEFAULT)
    snaps = list(range(t_start, t_last + 1, stride))
    file_out = FILE_OUT or ('./rrex_bottom_w_std_map.nc' if abs(dt - 43200.) < 1 else './rrex_bottom_w_std_map_%gh.nc' % (dt / 3600.))
    print('simulation %s: output every %g h, snapshots %d..%d every %d -> %d snapshots, written to %s'
          % (SIMUL, dt / 3600., t_start, t_last, stride, len(snaps), file_out))

    shape = np.asarray(sim.topo).shape
    keys = [None] + list(HAB_EXTRA)                           # None: lowest rho-level
    st_w, st_wu = {k: Welford(shape) for k in keys}, {k: Welford(shape) for k in keys}
    for n, t in enumerate(snaps):
        sim.update(time=t)
        if not has_snapshot(sim):
            print('missing or incomplete file for snapshot %d: %s' % (t, sim.ncfile)); continue
        for k, (wk, wuk) in bottom_fields(sim, HAB_EXTRA).items():
            st_w[k].add(wk); st_wu[k].add(wuk)
        if n % 25 == 0: print('snapshot %d (%d/%d)' % (t, n + 1, len(snaps)), flush=True)

    with Dataset(file_out, 'w') as nc:
        nc.createDimension('xi', shape[0]); nc.createDimension('eta', shape[1])
        nc.simulation, nc.dt_s, nc.stride = SIMUL, dt, stride
        nc.first_snapshot, nc.last_snapshot = t_start, t_last
        def put(name, arr, units=''):
            v = nc.createVariable(name, 'f4', ('xi', 'eta'), zlib=True); v[:] = np.asarray(arr, float); v.units = units
        nc.hab_extra = np.array(HAB_EXTRA, float)
        for k in keys:
            suf = '' if k is None else '_hab%g' % k
            put('std_w' + suf, st_w[k].std(), 'm day-1'); put('mean_w' + suf, st_w[k].mean_(), 'm day-1')
            put('std_wu' + suf, st_wu[k].std(), 'm day-1'); put('mean_wu' + suf, st_wu[k].mean_(), 'm day-1')
        put('count', st_w[None].n, '')
        put('h', sim.topo, 'm'); put('pm', sim.pm, 'm-1'); put('pn', sim.pn, 'm-1'); put('mask', sim.mask, '')
        put('lon', sim.x, 'degrees_east'); put('lat', sim.y, 'degrees_north')
    print('saved', file_out, '| snapshots used (max over the grid): %d/%d' % (st_w[None].n.max(), len(snaps)))


if __name__ == '__main__':
    main()
