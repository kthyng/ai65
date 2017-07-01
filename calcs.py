'''
Calculations.
'''

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmocean.cm as cmo
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import netCDF4 as netCDF
import argparse
from glob import glob
import numpy as np
from cartopy.io import shapereader
import matplotlib as mpl
mpl.rcParams.update({'font.size': 10})
import scipy.io
import os
import matplotlib.dates as Dates
# import octant
# import octant.depths
from plots import re

kind = 'baroclinic' #'barotropic'  # 'baroclinic'

locgrid = '/pong/raid/kthyng/froude/ai65/grid.nc'
locmodel = '/pong/raid/kthyng/froude/ai65/OUT/'
grid = netCDF.Dataset(locgrid)
pm = grid['pm'][0,0]; pn = grid['pn'][0,0]
dx = pm**-1; dy = pn**-1
h = grid['h'][:]
lon_psi = grid['lon_psi'][:]; lat_psi = grid['lat_psi'][:]
lon_rho = grid['lon_rho'][:]; lat_rho = grid['lat_rho'][:]
dates = pd.date_range("2006-09-01", "2006-10-01", freq="15T")
ntimes = len(dates)
g = 9.81  # m/s^2
Cd = 3e-3  # friction coefficient
rho0 = 1023.7

PEwest = np.empty(ntimes); KEwest = np.empty(ntimes)
PEnorth = np.empty(ntimes); KEnorth = np.empty(ntimes)
PEsouth = np.empty(ntimes); KEsouth = np.empty(ntimes)
Pin = np.empty(ntimes); Pout = np.empty(ntimes)
Pwest = np.empty(ntimes); Pnorth = np.empty(ntimes)
Psouth = np.empty(ntimes)
Pmix = np.empty(ntimes); Pfriction = np.empty(ntimes)
Pmom = np.empty(ntimes)
for i in range(ntimes):

    m = netCDF.Dataset(locmodel + 'ocean_his_' + str(i+1).zfill(4) + '.nc')

    ## Calculate over power dissipation across inlet in time
    # west
    ubarwest = re(m['ubar'][0,296:442,5:7],1).squeeze()  # interpolated to rho grid to match
    if kind == 'barotropic':
        vbarwest = re(m['vbar'][0,296:443,5],0).squeeze()
    elif kind == 'baroclinic':
        uwest = re(m['u'][0,:, 296:442,5:7],2).squeeze()
        vwest = re(m['v'][0,:, 296:443,5],1).squeeze()
    zetawest = m['zeta'][0,296:442,5]
    Hwest = (h[296:442,5] + zetawest)  # height
    PEwest[i] = rho0 * g * dy * (Hwest * ubarwest * zetawest).sum()  # potential energy, sum in y
    # PE = 0.5* rho0 * g * dy * (ubarwest * zetawest**2).sum()  # potential energy
    if kind == 'barotropic':
        KEwest[i] = 0.5 * rho0 * dy * (Hwest * ubarwest * (ubarwest**2 + vbarwest**2)).sum()  # sum in y
    elif kind == 'baroclinic':
        dzwest = Hwest/20.
        KEwest[i] = 0.5 * rho0 * dy * (dzwest * uwest * (uwest**2 + vwest**2)).sum()  # sum in z and y
    Pwest[i] = PEwest[i] + KEwest[i]

    # north
    vbarnorth = re(m['vbar'][0,440:442,5:38],0).squeeze()
    if kind == 'barotropic':
        ubarnorth = re(m['ubar'][0,440,5:39],0).squeeze()  # 0 here bc just vector
    elif kind == 'baroclinic':
        unorth = re(m['u'][0,:, 440,5:39],1).squeeze()
        vnorth = re(m['v'][0,:, 440:442,5:38],1).squeeze()
    zetanorth = m['zeta'][0,440,5:38]
    Hnorth = (h[440,5:38] + zetanorth)  # height
    PEnorth[i] = rho0 * g * dx * (Hnorth * vbarnorth * zetanorth).sum()  # potential energy
    # PE = 0.5 * rho0 * g * dx * (vbarnorth * zetanorth**2).sum()  # potential energy
    if kind == 'barotropic':
        KEnorth[i] = 0.5 * rho0 * dx  * (Hnorth * vbarnorth * (ubarnorth**2 + vbarnorth**2)).sum()  # sum in y
    elif kind == 'baroclinic':
        dznorth = Hnorth/20.
        KEnorth[i] = 0.5 * rho0 * dx * (dznorth * vnorth * (unorth**2 + vnorth**2)).sum()  # sum in z and y
    Pnorth[i] = PEnorth[i] + KEnorth[i]

    # east is closed

    # south
    vbarsouth = re(m['vbar'][0,5:7,130:-2],0).squeeze()
    if kind == 'barotropic':
        ubarsouth = re(m['ubar'][0,5,130:],0).squeeze()  # 0 here bc just vector
    elif kind == 'baroclinic':
        usouth = re(m['u'][0,:, 5,130:],1).squeeze()
        vsouth = re(m['v'][0,:, 5:7,130:-2],1).squeeze()
    zetasouth = m['zeta'][0,5,130:-2]
    Hsouth = (h[5,130:-2] + zetasouth)  # height
    PEsouth[i] = rho0 * g * dx * (Hsouth * vbarsouth * zetasouth).sum()  # potential energy
    # PE = 0.5 * rho0 * g * dx * (vbarsouth * zetasouth**2).sum()  # potential energy
    if kind == 'barotropic':
        KEsouth[i] = 0.5 * rho0 * dx * (Hsouth * vbarsouth * (ubarsouth**2 + vbarsouth**2)).sum()  # sum in y
    elif kind == 'baroclinic':
        dzsouth = Hsouth/20.
        KEsouth[i] = 0.5 * rho0 * dx * (dzsouth * vsouth * (usouth**2 + vsouth**2)).sum()  # sum in z and y
    Psouth[i] = PEsouth[i] + KEsouth[i]

    # Overall, then
    Pin[i] = Pwest[i] - Pnorth[i]  # west: + is in, for north: - is into domain
    Pout[i] = Psouth[i]  # Need to add this due to orientation
    ##
    # import pdb; pdb.set_trace()
    # ## Calculate mixing across inlet in time (buoyancy production)
    # AKs = m['AKs'][0, 1:-1, 5:-4, 5:-4].squeeze()  # vertical viscosity for tracer
    # rho = m['rho'][0, :, 5:-4, 5:-4].squeeze()  # density
    # zeta = m['zeta'][0,5:-4,5:-4]
    # dz = (h[5:-4,5:-4] + zeta)/20.  # easy since layers are uniform
    # drhodz = (rho[:-1] - rho[1:])/dz
    # Pmix[i] = g * dy * dx * (AKs*drhodz*dz).sum()
    # ##
    #
    # ## Calculate mixing across inlet in time (shear production)
    # AKv = m['AKv'][0, 1:-1, 5:-4, 5:-4].squeeze()  # vertical viscosity for momentum
    # u = re(m['u'][0, :, 5:-4, 5:-2],2)  # to rho grid
    # v = re(m['v'][0, :, 5:-2, 5:-4],1)  # to rho grid
    # dudz = (u[:-1] - u[1:])/dz
    # dvdz = (v[:-1] - v[1:])/dz
    # Pmom[i] = rho0 * dy * dx * (AKv*(dudz**2 + dvdz**2)*dz).sum()
    # ##
    #
    # ## Calculate power dissipation due to bottom friction, across domain
    # # using velocity from bottom of water column
    # Pfriction[i] = Cd * rho0 * dy * dx * (abs((u[0]**2 + v[0]**2)**(3/2))).sum()
    ##

df = pd.DataFrame(index=dates)
df['PEwest'] = PEwest
df['KEwest'] = KEwest
df['PEnorth'] = PEnorth
df['KEnorth'] = KEnorth
df['PEsouth'] = PEsouth
df['KEsouth'] = KEsouth
df['Pwest'] = Pwest
df['Pnorth'] = Pnorth
df['Psouth'] = Psouth
df['Pin'] = Pin
df['Pout'] = Pout
df['P'] = Pin + Pout
# df['Pmix'] = Pmix
# df['Pmom'] = Pmom
# df['Pfriction'] = Pfriction
df.to_csv('savedoutput/power_' + kind + '.csv')
# pd.read_csv('savedoutput/power.csv', index_col=0, parse_dates=True)

np.savez('savedoutput/power_' + kind + '_temp.npz', Pwest=Pwest, Pnorth=Pnorth, Psouth=Psouth, dates=dates)
