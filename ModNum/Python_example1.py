#!/usr/bin/env python

'''
Read data in a ROMS netcdf file

Plot a 2d horizontal field (temperature)

Save the figure

'''

#################################################
# Load Modules
#################################################

#for netcdf files
from netCDF4 import Dataset

#for numeric functions
import numpy as np

#plotting functions
import matplotlib.pyplot as py

#################################################
# Load data from netcdf file
#################################################

ncfile='his.nc'

#open netcdf file
nc = Dataset(ncfile, 'a')

# print all variables in file
print nc.variables

# print all dimensions
print nc.dimensions

#Load a 4d variable (default coordinates order= t,z,y,x)
rho=np.array(nc.variables['temp'])

# print variable size:
print rho.shape

#close netcdf file
nc.close()

#################################################
# Plot  data
#################################################

# Create a contour plot with 100 levels
py.contourf(rho[0,:,1,:],100,cmap=py.cm.spectral); 

# plot colorbar
py.colorbar()

# Save it as a .png file
py.savefig( 'image.png',  magnification='auto',bbox_inches='tight', dpi=300); 

# Show figure
py.show()
