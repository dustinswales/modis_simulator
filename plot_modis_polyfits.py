import os,numpy,netCDF4
import matplotlib.pyplot as plt

fileIN = 'modis_polyfits.nc'

dataIN      = netCDF4.Dataset(fileIN,'r')
re_liq      = dataIN.variables['re_liq'][:]
g_liq_new   = dataIN.variables['g_liq_new'][:]
g_liq_old   = dataIN.variables['g_liq_old'][:]
ssa_liq_new = dataIN.variables['ssa_liq_new'][:]
ssa_liq_old = dataIN.variables['ssa_liq_old'][:]
re_ice      = dataIN.variables['re_ice'][:]
g_ice_new   = dataIN.variables['g_ice_new'][:]
g_ice_old   = dataIN.variables['g_ice_old'][:]
ssa_ice_new = dataIN.variables['ssa_ice_new'][:]
ssa_ice_old = dataIN.variables['ssa_ice_old'][:]


# Plot optical-fields
fig, axs = plt.subplots(2, 2, figsize=(11, 7))
axs[0, 0].plot(re_liq,g_liq_new,'-')
axs[0, 0].plot(re_liq,g_liq_old,'--')
axs[0, 0].set_title('Liquid')
axs[0, 0].set_ylabel('Asymmetry parameter')
axs[0,0].legend(['new Coefficients','old coefficients'],loc='lower right')
axs[1, 0].plot(re_liq,ssa_liq_new,'-')
axs[1, 0].plot(re_liq,ssa_liq_old,'--')
axs[1, 0].set_ylabel('Single-scattering albedo')
axs[1, 0].set_xlabel('Re (microns)')
axs[0, 1].plot(re_ice,g_ice_new,'-')
axs[0, 1].plot(re_ice,g_ice_old,'--')
axs[0, 1].set_title('Ice')
axs[0, 1].set_ylabel('Asymmetry parameter')
axs[1, 1].plot(re_ice,ssa_ice_new,'-')
axs[1, 1].plot(re_ice,ssa_ice_old,'--')
axs[1, 1].set_ylabel('Single-scattering albedo')
axs[1, 1].set_xlabel('Re (microns)')
#plt.show()
plt.savefig('modis_polyfits.png')

