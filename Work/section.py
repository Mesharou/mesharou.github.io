#!/usr/bin/env python
###################################################################################

#import matplotlib
#matplotlib.use('Agg')

###################################################################################
# Load all useful modules (see /home/gula/Desktop/python_v2/Modules/)
###############################################
####################################

from Modules import *

import csv

from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText

###################################################################################
# Load simul
###################################################################################

#simul = load(simul = 'megatl6_n50_clim_newgrd_mean [700,1150,700,1450,[1,100,1]] 0')
#simul = load(simul = 'megatl6_n50_clim_newgrd_mean')
simul = load(simul = 'polgyr [750,1400,570,1750,[1,100,1]] 0')

fifig = './'

#####

save_netcdf = True
plot_temp = False

folder = './'

###################################################################################
# Read section of lat/lon
###################################################################################

#cr = csv.reader(open('OSNAP_EAST_line.csv','rb'))
cr = csv.reader(open('Angelina_North.csv','rb'))
cr = csv.reader(open('Angelina_South.csv','rb'))

lon = []; lat = []
for i,row in enumerate(cr): 
    if i>0:
        lon.append(float(row[0]))
        lat.append(float(row[1]))

###################################################################################
# find corresponding indices
###################################################################################


# gives only closest value
#xsec,ysec = tools.find_points(simul.x,simul.y, np.array(lon),np.array(lat))

[x,y] = np.mgrid[:simul.x.shape[0],:simul.x.shape[1]]

x0new  = interp.griddata((simul.x.ravel(),simul.y.ravel()),x.ravel(),(lon,lat),method='cubic')
y0new = interp.griddata((simul.x.ravel(),simul.y.ravel()),y.ravel(),(lon,lat),method='cubic')

###########################

imean = np.arange(np.floor(np.min(x0new))-1,np.ceil(np.max(x0new))+2)
poly = interp.interp1d(x0new,y0new,kind='linear',bounds_error=False,fill_value='extrapolate'); 

jmean = poly(imean)

#py.plot(x0new,y0new,'-o',); py.plot(imean,jmean,linewidth=3); py.show()

oldjmean = copy(jmean)

############################


x_0,sizex = np.int(imean[0]),-1

lonpsi = simul.x[x_0+1:sizex,:]
latpsi = simul.y[x_0+1:sizex,:]

################################################################################### 
# Define the new coordinates based on refined pos
###################################################################################

angle = sm.angle(jmean)

#Temporay test to cath the z=-50 isobath all the way
Ly = 10;

###########

py.plot(imean,jmean,'-g'); py.contourf(simul.topo.T,100); 
py.plot(imean,oldjmean,'-k',lw=4.); 
#py.plot(imean,jmean-Ly,'--g'); py.plot(imean,jmean+Ly,'--g'); 
py.savefig(fifig+'test.png',magnification='auto', dpi=250,bbox_inches='tight',transparence='true'); py.clf()

###########

[elem,coef,newXrot,newYrot] = sm.transform_grid(jmean,angle,simul.x[x_0+1:sizex,:].shape[1],Ly)

print 'transformation OK'

###################################################################################

newlon = np.sum(coef*lonpsi.ravel()[elem],2)
newlat = np.sum(coef*latpsi.ravel()[elem],2)

###################################################################################
#Plot transport
###################################################################################

toporot = np.sum(coef*simul.topo[x_0+1:sizex,:].ravel()[elem],2) #*maskrot
iy=Ly;

nx = toporot[:,iy].shape[0]

#####################

myvar = var('temp',simul).data

xstream = np.zeros((toporot.shape[0],myvar.shape[-1]+1))
xstream[:,-1] = newlon[:,iy]

for iz in range(myvar.shape[-1]):
    xstream[:,iz] = newlon[:,iy]

#####################



###################################################################################
# TIME LOOP
###################################################################################

for time in range(0,10000):
    simul.update(time)

    ##########

    [z_r,z_w] = tools.get_depths(simul)

    myvar = var('temp',simul).data
    salt = var('salt',simul).data
    u0 = tools.u2rho(var('u',simul).data)
    v0 = tools.v2rho(var('v',simul).data)

    for iz in range(u0.shape[-1]):
        [u0[:,:,iz],v0[:,:,iz]] = tools.rotuv(simul,u0[:,:,iz],v0[:,:,iz],psi=False)



    ###########

    temprot = np.zeros((toporot.shape[0],toporot.shape[1],myvar.shape[-1]))
    saltrot = np.zeros((toporot.shape[0],toporot.shape[1],myvar.shape[-1]))
    urot = np.zeros((toporot.shape[0],toporot.shape[1],myvar.shape[-1]))
    vrot = np.zeros((toporot.shape[0],toporot.shape[1],myvar.shape[-1]))
    zrot =  np.zeros((toporot.shape[0],toporot.shape[1],myvar.shape[-1]))


    nz = temprot.shape[-1]

    for iz in range(nz):
        temprot[:,:,iz] = np.sum(coef*myvar[x_0+1:sizex,:,iz].ravel()[elem],2)
        saltrot[:,:,iz] = np.sum(coef*salt[x_0+1:sizex,:,iz].ravel()[elem],2)
        urot[:,:,iz] = np.sum(coef*u0[x_0+1:sizex,:,iz].ravel()[elem],2)
        vrot[:,:,iz] = np.sum(coef*v0[x_0+1:sizex,:,iz].ravel()[elem],2)
        zrot[:,:,iz] = np.sum(coef*z_r[x_0+1:sizex,:,iz].ravel()[elem],2)



    if save_netcdf:
        ###################################################################################
        # save netcdf file
        ###################################################################################

        filename = folder + simul.simul + '_OSNAP.' + '{0:04}'.format(simul.filetime) + '.nc'



        if simul.infiletime==0:

            nc = Dataset(filename,'w',format='NETCDF4')
            # --- create global attributes ------------------------------------------
            nc.title        = 'OSNAP Section'

            nc.createDimension('time',None)
            nc.createDimension('lon',nx)
            nc.createDimension('sig',nz)

            netvar             = nc.createVariable('lon','d',('lon',))
            netvar.long_name   = 'longitude'
            netvar             = nc.createVariable('sig','d',('sig',))
            netvar.long_name   = 'sigma level'

            nc.variables['lon'][:] = newlon[:,iy]
            nc.variables['sig'][:] = range(nz)

            netvar             = nc.createVariable('temp','d',('time','sig','lon',))
            netvar.long_name   = 'temperature'
            netvar.units       = 'deg'

            netvar             = nc.createVariable('salt','d',('time','sig','lon',))
            netvar.long_name   = 'salinity'
            netvar.units       = 'PSU'

            netvar             = nc.createVariable('u','d',('time','sig','lon',))
            netvar.long_name   = 'zonal velocity'
            netvar.units       = 'm/s'

            netvar             = nc.createVariable('v','d',('time','sig','lon',))
            netvar.long_name   = 'meridional velocity'
            netvar.units       = 'm/s'

            netvar             = nc.createVariable('z','d',('time','sig','lon',))
            netvar.long_name   = 'depth'
            netvar.units       = 'm'

        else:

            nc = Dataset(filename,'a',format='NETCDF4')

        nc.variables['temp'][simul.infiletime,:,:] = temprot[:,iy,:].T
        nc.variables['salt'][simul.infiletime,:,:] = saltrot[:,iy,:].T
        nc.variables['u'][simul.infiletime,:,:] = urot[:,iy,:].T
        nc.variables['v'][simul.infiletime,:,:] = vrot[:,iy,:].T
        nc.variables['z'][simul.infiletime,:,:] = zrot[:,iy,:].T

        nc.close()

    ###################################################################################

    if plot_temp:
        ###################################################################################
        # PLOT
        ###################################################################################

        my_cmap=py.cm.jet
        levels= np.arange(2.,10.,0.1)    
        varunit = 'T [$^{\circ}$C]'

        ###################################################################################
        # PLOT
        ###################################################################################

        cbarlabelsize = 6
        fontsize0 = 12
        fontsize1 = 16
        fontsize2 = format(22)
        fontsize3 = format(28)
        font = {'size'   : fontsize1}
        py.rc('font', **font)

        ################################

        fig = py.figure(figsize=(9.0,3.0))
        ax = fig.add_subplot(111,axisbg='Gainsboro')

        ################################

        joe = ax.contourf(xstream[:,:-1],zrot[:,iy,:],temprot[:,iy,:],levels,extend='both',cmap=my_cmap)

        #joe = ax.pcolormesh(xstream[:,:-1],zrot[:,iy,:],ma.masked_invalid(temprot[:,iy,:]),vmin=levels.min(),vmax=levels.max(),cmap = my_cmap,rasterized=True); 

        for i in range(len(lon)):
            py.plot([lon[i],lon[i]],[-3050,-3000],'-k',lw=2.)

        ax.plot(xstream[:,0],-toporot[:,iy],'-k',lw=3.)
        py.ylabel('z [m]',fontsize=fontsize2)

        #ax.set_ylim([-1500., 0.])
        #py.xlim([11., 59.])

        #######################

        at = AnchoredText(varunit,loc=4, prop=dict(size=fontsize2), frameon=True, )
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        ax.add_artist(at)

        py.subplots_adjust(right=0.94)
        cbar_ax = fig.add_axes([0.96, 0.2, 0.02, 0.60])
        fig.colorbar(joe, cax=cbar_ax)

        ########################################################################

        py.savefig(fifig+ 'vertical_section_temp_300.png',magnification='auto', dpi=300,bbox_inches='tight',transparence='true');

        ###################################################################################


































































