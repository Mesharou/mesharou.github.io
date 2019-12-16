# ----------------------------------------------------------------------
# Rotating Shallow Water model
#
# author: G. Roullet
# email: roullet@univ-brest.fr
# date: november 2018
#
# Modified by J. Gula / January 2019
#
# The model uses a C-grid in space
# It is formulated in vector-invariant form
#
# The domain is closed with no flux+free-slip boundary condition
#
# Several discretizations are possible form
# - Coriolis term: centered, upwind 1st, 3rd or 5th order
# - mass flux: centered, upwind 1st
#
# Upwinded discretizations introduce a numerical dissipation that
# helps having smoother solutions, especially on the PV
#
# Several time schemes are proposed
#
# Three initial conditions are proposed: a dam-break problem,
# a geostrophic adjustment and a vortex merging
#
# The RHS computation needs the 'jit' decorator from the numba module
# you need to have numba installed.
#
# Results are shown on the fly. You may plot 'pv', 'u', 'v' or 'h'
# If you have 'avconv' or 'ffmeg' installed, you may set generate_mp4 = True
# to have mp4 generated on the fly
#

import matplotlib
from matplotlib.ticker import NullFormatter

macuser = False
matplotlib.use('TkAgg')
font = {'size': 16}
matplotlib.rc('font', **font)

import subprocess
import os
from numba import jit
import numpy as np
import matplotlib.pyplot as plt

plt.ion()

Lx, Ly = 10, 10  # domain lengths
nx, ny = 100, 100  # resolution (in grid cells)

g = 1.  # acceleration of gravity
Hmax = 1.  # total depth at rest (prognostic 'hC' is the depth anomaly)
f = 0.  # Coriolis parameter
#tend = 1. # integration time
tend = 2. * Lx / np.sqrt(g* Hmax) # integration time (defined as time for a gravity wave to travel the size of the domain)

print('tend is ' + str(tend))
cfl = 1.587/2
# cfl = 1.75/2
# cfl = 0.4

timescheme = 'RK3_SSP'  #
vort_discretiz = 'upwind'
vort_order = 5  # order is 1, 3 or 5, matters only when vort_discretiz = 'upwind'
cont_discretiz = 'centered'  # can be 'centered' or 'upwind'

umax = np.sqrt(g*Hmax)
dx, dy = Lx/nx, Ly/ny
dt = cfl/np.sqrt((umax/dx)**2+(umax/dy)**2)
nt = int(tend/dt)

# plotting options
plotperiod = 1  # refresh plot every plotperiod iterations
generate_mp4 = True  # turn it off if you don't have 'avconv' or 'ffmeg' installed
cax = [-0.0005, 0.0005]  # min-max values for the colorbar
colorscheme = 'auto'  # if 'imposed' use 'cax', if 'auto' or 'auto_sym' to self adjust
plotvar = 'h'


# Type of topography: ['flat', 'jump', etc.]
topography = 'inlet';

# Initial condition ['geostrophicadj',etc.]
initial_condition = 'geostrophicadj'

# Figure/Movie name
case = initial_condition + '_' + topography + '_' + 'H01_forcing'

#####################

uU = np.zeros((ny, nx+1))
vV = np.zeros((ny+1, nx))
hC = np.zeros((ny, nx))
H = np.zeros((ny, nx))

# coordinates of the cells centers 'C' and cell corners 'F'
xC = (0.5+np.arange(nx))*dx
xF = np.arange(nx+1)*dx
yC = (0.5+np.arange(ny))*dy
yF = np.arange(ny+1)*dy

def ddx(z):
    return (z[:, 1:]-z[:, :-1])/dx

def ddy(z):
    return (z[1:, :]-z[:-1, :])/dy

def avx(z):
    return (z[:, 1:]+z[:, :-1])*0.5

def avy(z):
    return (z[1:, :]+z[:-1, :])*0.5

def curlF(U_, V_):
    omega = np.zeros((ny+1, nx+1))
    omega[1:-1, :] = ddy(U_)
    omega[:, 1:-1] -= ddx(V_)
    return omega


def curlC(U_, V_):
    return avx(avy(curlF(U_, V_)))


def pvC(h_, U_, V_):
    return (f+curlC(U_, V_))/(H+h_)


@jit
def rhs(state):
    """ compute the RHS of the RSW equations
    input: state is a list of three arrays
    output: a list of three arrays with the tendencies

    The code is written with explicit i,j loops. @jit compiles
    this routine prior the execution. It makes the code as fast as if
    it where coded directly without loops"""
    hC_, uU_, vV_ = state
    
    
    

    if cont_discretiz == 'upwind':
        dhC = np.zeros_like(hC_)
        # the mass flux is upwinded (1st order)
        cff = 1/dx
        for j in range(ny):
            for i in range(1, nx):
                fu = uU_[j, i]
                if fu > 0:
                    flux = (H[j, i-1]+hC_[j, i-1])*fu*cff
                else:
                    flux = (H[j, i]+hC_[j, i])*fu*cff
                dhC[j, i-1] -= flux
                dhC[j, i] += flux

        cff = 1/dy
        for j in range(1, ny):
            for i in range(nx):
                fv = vV_[j, i]
                if fv > 0:
                    flux = (H[j-1, i]+hC_[j-1, i])*fv*cff
                else:
                    flux = (H[j, i]+hC_[j, i])*fv*cff
                dhC[j-1, i] -= flux
                dhC[j, i] += flux

    elif cont_discretiz == 'centered':
        # if jit is not available, use this centered implementation
        # centered mass flux
        hU = np.zeros_like(uU_)
        hV = np.zeros_like(vV_)
        hU[:, 1:-1] = avx(H+hC_)
        hV[1:-1, :] = avy(H+hC_)

        dhC = -(ddx(hU*uU_) + ddy(hV*vV_))

    else:
        raise ValueError('cont_discretiz should be "centered" or "upwind"')

    # vorticity at corner cell
    omega = f+curlF(uU_, vV_)

    # Bernoulli function (square then average, not the reverse)
    bC = g*hC_+0.5*(avx(uU_**2) + avy(vV_**2))
    # bC = g*hC_+0.5*(avx(uU_)**2 + avy(vV_)**2)

    duU = np.zeros_like(uU_)
    dvV = np.zeros_like(vV_)

    # coefficients for 3rd order biased interpolation
    c1, c2, c3 = -1/6., 5/6., 2/6.
    # and for the 5th order (volume flux coefficients)
    d1, d2, d3, d4, d5 = 1./30., -13./60., 47./60., 9./20., -1./20.

    if vort_discretiz == 'upwind':
        FV = avy(avx(vV_))
        FU = avx(avy(uU_))
        for j in range(ny):
            for i in range(nx-1):
                if (vort_order == 5) and (j < ny-2) and (j > 0):
                    op = (+ d5*omega[j-1, i+1]
                          + d4*omega[j, i+1]
                          + d3*omega[j+1, i+1]
                          + d2*omega[j+2, i+1]
                          + d1*omega[j+3, i+1])
                elif (vort_order >= 3) and (j < ny-1):
                    op = (+ c3*omega[j, i+1]
                          + c2*omega[j+1, i+1]
                          + c1*omega[j+2, i+1])
                else:
                    op = omega[j+1, i+1]

                if (vort_order == 5) and (j < ny-1) and (j > 1):
                    om = (+ d1*omega[j-2, i+1]
                          + d2*omega[j-1, i+1]
                          + d3*omega[j, i+1]
                          + d4*omega[j+1, i+1]
                          + d5*omega[j+2, i+1])
                if (vort_order >= 3) and (j > 0):
                    om = (+ c1*omega[j-1, i+1]
                          + c2*omega[j, i+1]
                          + c3*omega[j+1, i+1])
                else:
                    om = omega[j, i+1]
                FV[j, i] = max(FV[j, i], 0)*op + min(FV[j, i], 0)*om

        for j in range(ny-1):
            for i in range(nx):
                if (vort_order == 5) and (i < nx-2) and (i > 0):
                    op = (+ d5*omega[j+1, i-1]
                          + d4*omega[j+1, i]
                          + d3*omega[j+1, i+1]
                          + d2*omega[j+1, i+2]
                          + d1*omega[j+1, i+3])
                elif (vort_order >= 3) and (i < nx-1):
                    op = (+ c3*omega[j+1, i]
                          + c2*omega[j+1, i+1]
                          + c1*omega[j+1, i+2])
                else:
                    op = omega[j+1, i+1]

                if (vort_order == 5) and (i < nx-1) and (i > 1):
                    om = (+ d1*omega[j+1, i-2]
                          + d2*omega[j+1, i-1]
                          + d3*omega[j+1, i]
                          + d4*omega[j+1, i+1]
                          + d5*omega[j+1, i+2])
                elif (vort_order >= 3) and (i > 0):
                    om = (+ c1*omega[j+1, i-1]
                          + c2*omega[j+1, i]
                          + c3*omega[j+1, i+1])
                else:
                    om = omega[j+1, i]
                FU[j, i] = max(FU[j, i], 0)*op + min(FU[j, i], 0)*om

        # if jit is not available, use this 1st order upwind implementation
        # FV = (FV-np.abs(FV))*omega[:-1, 1:-1] + (FV+np.abs(FV))*omega[1:, 1:-1]
        # FU = (FU-np.abs(FU))*omega[1:-1, :-1] + (FU+np.abs(FU))*omega[1:-1, 1:]
        # FV *= 0.5
        # FU *= 0.5

        duU[:, 1:-1] = -ddx(bC) + FV
        dvV[1:-1, :] = -ddy(bC) - FU

    elif vort_discretiz == 'centered':
        # centered discretization
        duU[:, 1:-1] = -ddx(bC) + avy(omega[:, 1:-1])*avy(avx(vV_))
        dvV[1:-1, :] = -ddy(bC) - avx(omega[1:-1, :])*avx(avy(uU_))

        # another one
        # duU[:, 1:-1] = -ddx(bC) + avy(omega[:, 1:-1])*avy(avx(vV_))
        # dvV[1:-1, :] = -ddy(bC) - avx(omega[1:-1, :])*avx(avy(uU_))

    else:
        raise ValueError('vort_discretiz should be "centered" or "upwind"')

    return [dhC, duU, dvV]






def dambreak(x1, y1):
    xx, yy = np.meshgrid(x1, y1)
    # set a non zero slope to have a slanted dam break
    slope = 0.
    sigma = 0.08
    return np.tanh(((xx-Lx/2)-slope*(yy-Ly/2))/sigma)

def vortex(x1, y1, x0, y0, d):
    xx, yy = np.meshgrid(x1, y1)
    d2 = (xx-x0)**2 + (yy-y0)**2
    return np.exp(-d2/(2*d**2))

def jump(x1, y1, x0, y0, d):
    xx, yy = np.meshgrid(x1, y1)
    h = np.zeros_like(xx)
    h[np.logical_and(xx>x0-d,xx<x0+d)] = 1.
    return h

###############################
# Topography
###############################

if topography == 'flat':
    H[:, :] = Hmax
elif topography == 'jump':
    H[:, :nx//2] = Hmax
    H[:, nx//2:] = Hmax/10.
elif topography == 'inlet':
    H[:, 2*nx//3:] = 1.
    H[:ny//2-10, :] = 1e-5
    H[ny//2+10:, :] = 1e-5
    H[:, :nx//3] = Hmax
    H[:, nx//3:2*nx//3] = ( Hmax * (xC[nx//3:2*nx//3] - xC[nx//3]) + (xC[2*nx//3] - xC[nx//3:2*nx//3]) ) / (xC[2*nx//3] - xC[nx//3])
elif topography == 'slope':
    xx, yy = np.meshgrid(xC, yC)
    H[:, :] = Hmax * (1. - xC/xC.max())
elif topography == 'bowl':
    H[:, :] = Hmax * vortex(xC, yC, 0.5*Lx, 0.5*Ly, Lx/5)
else:
    H[:, :] = Hmax


plt.imshow(H); plt.show(); plt.pause(5)

###############################
# Initial Condition
###############################


rossby = 0.

if initial_condition == 'dambreak':
    # use a small deformation radius sqrt(g*H)/f
    # 0.1 of Lx to clearly see the Kelvin waves
    # propagating along the boundary
    hC0 = 0.1*dambreak(xC, yC)

if initial_condition == 'geostrophicadj':
    d = dx
    amp = 0.01
    hC0 = amp*vortex(xC, yC, 0.5*Lx, 0.5*Ly, d)


if initial_condition == 'vortexmerging':
    d = 0.1
    # the vortex amplitude controls the Froude number
    amp = 0.05
    hC0 = amp*vortex(xC, yC, 0.5, 0.4, d)
    hC0[:, :] += amp*vortex(xC, yC, 0.5, 0.6, d)

    # to set initial geostropic adjustement
    # define exactly the same height but at corner cells...
    hF = amp*vortex(xF, yF, 0.5, 0.4, d)
    hF[:, :] += amp*vortex(xF, yF, 0.5, 0.6, d)
    # then take the rotated gradient of it
    uU = -(g/f)*ddy(hF)
    vV = +(g/f)*ddx(hF)

    maxspeed = max(np.max(uU.ravel()), np.max(vV.ravel()))
    froude = maxspeed/np.sqrt(g*Hmax)
    rossby = maxspeed/(f*d)
    if f>0:
        print("Froude = %.2f" % froude)
        print("Rossby = %.2f" % rossby)


class Movie(object):
    """ Home made class to generate mp4 """

    def __init__(self, fig, name=case + 'mymovie'):
        """ input: fig is the handle to the figure """
        self.fig = fig
        canvas_width, canvas_height = self.fig.canvas.get_width_height()
        # Open an ffmpeg process
        outf = '%s.mp4' % name
        videoencoder = None
        for v in ['avconv', 'ffmeg']:
            if subprocess.call(['which', v], stdout=subprocess.PIPE) == 0:
                videoencoder = v

        if videoencoder is None:
            print('\n')
            print('Neither avconv or ffmeg was found')
            print('Install one or set param.generate_mp4 = False')
            raise ValueError('Install avconv or ffmeg')

        cmdstring = (videoencoder,
                     '-y', '-r', '30',  # overwrite, 30fps
                     # size of image string
                     '-s', '%dx%d' % (canvas_width, canvas_height),
                     '-pix_fmt', 'argb',  # format
                     '-f', 'rawvideo',
                     # tell ffmpeg to expect raw video from the pipe
                     '-i', '-',
                     '-vcodec', 'libx264', outf)  # output encoding

        devnull = open(os.devnull, 'wb')
        self.process = subprocess.Popen(cmdstring,
                                        stdin=subprocess.PIPE,
                                        stdout=devnull,
                                        stderr=devnull)

    def addframe(self):
        string = self.fig.canvas.tostring_argb()
        self.process.stdin.write(string)

    def finalize(self):
        self.process.communicate()


class Figure(object):
    def __init__(self, varname):
        self.plotvar = varname
        z2d = getplotvar(self.plotvar)
        
        self.fig = plt.figure(figsize=(8, 8))
        #gs1 = self.fig.add_gridspec(nrows=5, ncols=1,  wspace=0.01)

        #self.ax1 = self.fig.add_subplot(gs1[:-1, :])
        #self.ax1 = self.fig.add_subplot(1, 1, 1)
        self.ax1 = plt.subplot2grid((5,1),(0,0),rowspan=4)  

        self.im = self.ax1.imshow(z2d, extent=[0, Lx, 0, Ly], cmap='RdBu_r',interpolation='nearest', origin='lower')
        if colorscheme == 'imposed':
            plt.colorbar(self.im,shrink=0.5)
        self.ax1.set_title('%s / t=%.2f' % (varname, time))
        self.ax1.set_ylabel('Y')
        if f>0:
            self.ax1.text(0.02*Lx, 0.95*Ly, 'Rd=%.2f' % (np.sqrt(g*Hmax)/f))
            self.ax1.text(0.02*Lx, 0.9*Ly, 'Ro=%.2f' % rossby)
        self.ax1.text(0.02*Lx, 0.85*Ly, r'$\sqrt{gH}$=%.2f' % (np.sqrt(g*Hmax)))
        
        #self.ax2 = self.fig.add_subplot(gs1[-1, :])
        #self.ax2 = self.fig.add_subplot(10, 1, 10)
        self.ax2 = plt.subplot2grid((6,1),(-1,0))  

        self.im2, = self.ax2.plot(xC, z2d[ny//2,:]);
        self.ax2.yaxis.set_major_formatter( NullFormatter() )
        self.ax2.set_xlabel('X')
        self.fig.show()
        if macuser:
            plt.pause(0.0001)
        else:
            self.fig.canvas.draw()
        if generate_mp4:
            self.mov = Movie(self.fig)
        self.update()

    def update(self):
        z2d = getplotvar(self.plotvar)
        self.im.set_array(z2d)
        if colorscheme == 'imposed':
            vmin, vmax = cax
        elif colorscheme == 'auto':
            vmin, vmax = np.min(z2d), np.max(z2d)
        elif colorscheme == 'auto_sym':
            vmax = np.max(np.abs(z2d))
            vmin= -vmax
        else:
            raise ValueError('colorscheme should be "imposed" or "auto"')

        self.ax1.set_title('%s / t=%.2f' % (self.plotvar, time))
        self.im.set_clim(vmin=vmin, vmax=vmax)
        
        self.im2.set_ydata(z2d[ny//2,:])
        if 'auto' in colorscheme:
            self.ax2.set_ylim([vmin,vmax])
        
        if macuser:
            plt.pause(0.0001)
        else:
            self.fig.canvas.draw()
        if generate_mp4:
            self.mov.addframe()


def getplotvar(plotvar):
    """ return the desired 2D array

    Assumes that pv, uU, vV, hC are global variables
    But, unlike a classic Fortran code, these arrays
    are created/destroyed continuously in the time loop"""
    if plotvar == 'pv':
        z2d = pv
    elif plotvar == 'u':
        z2d = uU
    elif plotvar == 'v':
        z2d = vV
    elif plotvar == 'h':
        z2d = hC
    else:
        raise ValueError('plotvar should be "pv", "u", "v" or "h"')
    return z2d


Volume = np.zeros((nt,))
Energy = np.zeros((nt,))
Potene = np.zeros((nt,))
PVcont = np.zeros((nt,))

hC[:,:] = hC0

state = [hC, uU, vV]
pv = pvC(hC, uU, vV)


time = 0.
plt.close('all')
fig = Figure(plotvar)

np.savez(case + 'initialstate', arr1=xC,arr2=yC,arr3=hC)

for kt in range(nt):
    #print('write = %i/%i -- time = %5.2f/%.2f' % (kt, nt, time, tend), end='')
    
    if timescheme == 'Heun':
        ds0 = rhs(state)
        s0 = [x+y*dt for x, y in zip(state, ds0)]

        ds1 = rhs(s0)
        state = [x+(y+z)*(dt/2) for x, y, z in zip(state, ds0, ds1)]

    if timescheme == 'LeapFrog':
        ds0 = rhs(state)
        if kt == 0:
            sb = state
            state = [x+y*dt for x, y in zip(state, ds0)]
        else:
            s1 = [x+y*(2*dt) for x, y in zip(sb, ds0)]
            sb = state
            state = s1

    if timescheme == 'RK3_SSP':
        # SSP RK3
        ds0 = rhs(state)
        s0 = [x+y*dt for x, y in zip(state, ds0)]

        ds1 = rhs(s0)
        s1 = [x+(y+z)*(dt/4) for x, y, z in zip(state, ds0, ds1)]

        ds2 = rhs(s1)
        state = [x+(w+y+4*z)*(dt/6)
                 for x, w, y, z in zip(state, ds0, ds1, ds2)]

    if timescheme == 'LFAM3':
        # this is ROMS/CROCO time scheme
        # it's third order and seems to compare quite well with RK3_SSP
        ds0 = rhs(state)
        if kt == 0:
            sb = state
            state = [x+y*dt for x, y in zip(state, ds0)]
        else:
            s1 = [x+y*(2*dt) for x, y in zip(sb, ds0)]
            s2 = [(5*x + 8*y-z)/12 for x, y, z in zip(s1, state, sb)]
            ds1 = rhs(s2)
            sb = state
            state = [x+y*dt for x, y in zip(state, ds1)]

    hC, uU, vV = state
    pv[:, :] = pvC(hC, uU, vV)

    time += dt
    energy = 0.5*(g*(2*H*hC+hC**2)+(H+hC)*(avx(uU**2)+avy(vV**2)))
    pe = 0.5*(g*(2*H*hC+hC**2))

    if kt % plotperiod == 0:
        fig.update()

    Volume[kt] = np.mean(hC.ravel())
    Energy[kt] = np.mean(energy.ravel())
    Potene[kt] = np.mean(pe.ravel())
    PVcont[kt] = np.mean(pv.ravel())

print('/done')
fig.update()
if generate_mp4:
    fig.mov.finalize()
plt.savefig(case + 'finalstate.pdf')

np.savez(case + 'finalstate', arr1=xC,arr2=yC,arr3=hC)

plt.figure(figsize=(10, 6))
time = np.arange(nt)*dt
plt.plot(time, (Energy-Energy[0])/Energy[0], label=r'$\Delta E\,/\,E_0$')
plt.plot(time, (Potene-Potene[0])/Potene[0], label=r'$\Delta P\,/\,P_0$')
#plt.plot(time, (PVcont-PVcont[0])/PVcont[0], label=r'$\Delta q\,/\,q_0$')
plt.legend()
plt.xlabel('time')
plt.show()













