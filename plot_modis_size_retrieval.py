import os,numpy,netCDF4
import matplotlib.pyplot as plt
from matplotlib import cm

fileIN = 'modis_polyfits.nc'

dataIN      = netCDF4.Dataset(fileIN,'r')
re_liq      = dataIN.variables['re_liq'][:]
re_ice      = dataIN.variables['re_ice'][:]
tau         = dataIN.variables['tau'][:]
sizeLIQ_old = dataIN.variables['sizeLIQ_old'][:,:]
sizeLIQ_new = dataIN.variables['sizeLIQ_new'][:,:]
sizeICE_old = dataIN.variables['sizeICE_old'][:,:]
sizeICE_new = dataIN.variables['sizeICE_new'][:,:]

# Create masks
sizeLIQ_oldm = sizeLIQ_old > 0.
sizeLIQ_newm = sizeLIQ_new > 0.
sizeICE_oldm = sizeICE_old > 0.
sizeICE_newm = sizeICE_new > 0.

sizeLIQ_old = numpy.where(sizeLIQ_old < 0., 0., sizeLIQ_old)
sizeLIQ_new = numpy.where(sizeLIQ_new < 0., 0., sizeLIQ_new)
sizeICE_old = numpy.where(sizeICE_old < 0., 0., sizeICE_old)
sizeICE_new = numpy.where(sizeICE_new < 0., 0., sizeICE_new)

# Plot retrieval mask
fig, axs = plt.subplots(2, 2, figsize=(7, 7))
axs[0,0].imshow(sizeLIQ_oldm,extent=(min(tau),max(tau),max(re_liq),min(re_liq)), \
                cmap='binary', aspect='auto')
axs[0,0].set_title("Size-liquid (old)")
axs[0,0].set_ylabel("Size (microns)")
axs[0,1].imshow(sizeLIQ_newm,extent=(min(tau),max(tau),max(re_liq),min(re_liq)), \
                cmap='binary', aspect='auto')
axs[0,1].set_title("Size-liquid (new)")
axs[1,0].imshow(sizeICE_oldm,extent=(min(tau),max(tau),max(re_ice),min(re_ice)), \
                cmap='binary', aspect='auto')
axs[1,0].set_title("Size-ice (old)")
axs[1,0].set_xlabel("Optical-depth")
axs[1,0].set_ylabel("Size (microns)")
axs[1,1].imshow(sizeICE_newm,extent=(min(tau),max(tau),max(re_ice),min(re_ice)), \
                cmap='binary', aspect='auto')
axs[1,1].set_title("Size-ice (new)")
axs[1,1].set_xlabel("Optical-depth")
#plt.show()
plt.savefig('modis_size_retrieval.png')
