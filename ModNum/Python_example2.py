#!/usr/bin/env python



#################################################
# Load Modules
#################################################

#for netcdf files
from netCDF4 import Dataset

#for numeric functions
import numpy as np

#plotting functions
import matplotlib.pyplot as py
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter

#############

import sys
sys.path.append("./Modules/")

#some ROMS tools examples 
import R_tools as tools
import R_tools_fort as toolsF


###################################################################################
# Use fortran order for arrays
####################################################################################

def Forder(var):
   return np.asfortranarray(var.T,dtype=np.float64)

####################################################################################


# Choose a time-step

itime=10


#################################################
# Load data from netcdf file
#################################################


ncfile='roms_his.nc'

#open netcdf file
nc = Dataset(ncfile, 'a')


#Load variables
rho=Forder(nc.variables['rho'][itime,:,:,:])
u=Forder(nc.variables['u'][itime,:,:,:])
v=Forder(nc.variables['v'][itime,:,:,:])

#Load some parameters

zeta=Forder(nc.variables['zeta'][itime,:,:])
topo=Forder(nc.variables['h'][:])
pm=Forder(nc.variables['pm'][:])
pn=Forder(nc.variables['pn'][:])

hc = nc.hc
Cs_r = Forder(nc.Cs_r)
Cs_w = Forder(nc.Cs_w)

#close netcdf file
nc.close()

###################################################################################
#Compute vertical coordinates 
###################################################################################


(z_r,z_w) = toolsF.zlevs(topo, zeta[:,:], hc, Cs_r, Cs_w)
[x,y,z] =np.mgrid[0:z_r.shape[0],0:z_r.shape[1],0:z_r.shape[2]]
x = (x.T/pm.T).T * 1e-3

[x_w,y0,z0] =np.mgrid[0:z_w.shape[0],0:z_w.shape[1],0:z_w.shape[2]]
x_w = (x_w.T/pm.T).T * 1e-3

#################################################
#Compute vertical velocity
#################################################

w = toolsF.get_wvlcty(u,v,z_r,z_w,pm,pn)

#################################################
#Compute dwdz
#################################################

dwdz = (w[:,:,1:] - w[:,:,:-1])/(z_r[:,:,1:] - z_r[:,:,:-1])

#################################################
#Compute v * w dwdz
#################################################


vwdwdz = tools.v2rho(tools.rho2w(v)) * tools.rho2w(w) * dwdz


#################################################
# Plot  data
#################################################

fig = py.figure(figsize=(8.0,6.0)) 
fig.set_tight_layout(2.)

###################

ax1 = py.subplot(2,1,1);
py.contourf(x[:,1,:],z_r[:,1,:],w[:,1,:],100,cmap=py.cm.spectral);py.colorbar()
py.contour(x[:,1,:],z_r[:,1,:],rho[:,1,:],25,  linewidths = (1.,));
CS1 = py.contour(x[:,1,:],z_r[:,1,:],np.max(z) - z[:,1,:],range(100),colors = ('k',), linewidths = (0.2,)); 
py.clabel(CS1, fmt = '%2.0f', colors = 'k', fontsize=6)  ; 
ax1.xaxis.set_major_formatter( NullFormatter() )
py.ylabel(r'$z\,(m)$',fontsize=18)
py.title(r'$w$', fontsize=20)

##############

py.subplot(2,1,2);
py.contourf(x_w[:,1,1:-1],z_w[:,1,1:-1],vwdwdz[:,1,:],100,cmap=py.cm.spectral);py.colorbar()
py.contour(x[:,1,:],z_r[:,1,:],rho[:,1,:],25,  linewidths = (1.,));
CS1 = py.contour(x[:,1,:],z_r[:,1,:],np.max(z) - z[:,1,:],range(100),colors = ('k',), linewidths = (0.2,)); 
py.clabel(CS1, fmt = '%2.0f', colors = 'k', fontsize=6)  ; 
py.xlabel(r'$x\,(km)$',fontsize=18)
py.ylabel(r'$z\,(m)$',fontsize=18)
py.title(r'$v w \frac{dw}{dz}$',verticalalignment='bottom', fontsize=20)

# Save it as a .png file
py.savefig( 'image.png',  magnification='auto',bbox_inches='tight', dpi=300); 
#
# Show figure
py.show()

#################################################



