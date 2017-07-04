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
from plots import re

which = 'all' # 'all' 'hor'
kind = 'barotropic'
doterms = False

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


if which == 'hor':
    # west transects: [0:1, 2]
    twest = np.array([[245, 340, 65], [0, 245, 160]])
    # east transects: [0:1, 2]
    teast = np.array([[310, 340, 120], [0, 310, 200]])
    # north transects: [0, 1:2]
    tnorth = np.array([[340, 65, 120], [310,120,200]])
    # south transects: [0, 1:2]
    tsouth = np.array([[245, 65, 160]])
elif which == 'all':
    twest = np.array([[296, 442, 5]])
    teast = []
    tnorth = np.array([[440, 5, 38]])
    tsouth = np.array([[5, 130, 295]])

PEeast = np.empty(ntimes); KEeast = np.empty(ntimes)
PEwest = np.empty(ntimes); KEwest = np.empty(ntimes)
PEnorth = np.empty(ntimes); KEnorth = np.empty(ntimes)
PEsouth = np.empty(ntimes); KEsouth = np.empty(ntimes)
Pin = np.empty(ntimes); Pout = np.empty(ntimes)
Pwest = np.empty(ntimes); Pnorth = np.empty(ntimes)
Psouth = np.empty(ntimes); Peast = np.empty(ntimes)
if doterms:
    Pmix = np.empty(ntimes); Pfriction = np.empty(ntimes)
    Pmom = np.empty(ntimes)
for i in range(ntimes):

    m = netCDF.Dataset(locmodel + 'ocean_his_' + str(i+1).zfill(4) + '.nc')
    # Calculate over power dissipation across inlet in time
    for vec in twest:  # west
        ubarwest = re(m['ubar'][0,vec[0]:vec[1],vec[2]:vec[2]+2],1).squeeze()  # interpolated to rho grid to match
        if kind == 'barotropic':
            vbarwest = re(m['vbar'][0,vec[0]:vec[1]+1,vec[2]],0).squeeze()
        elif kind == 'baroclinic':
            uwest = re(m['u'][0,:, vec[0]:vec[1],vec[2]:vec[2]+2],2).squeeze()
            vwest = re(m['v'][0,:, vec[0]:vec[1]+1,vec[2]],1).squeeze()
        zetawest = m['zeta'][0,vec[0]:vec[1],vec[2]]
        Hwest = (h[vec[0]:vec[1],vec[2]] + zetawest)  # height
        # PEwest[i] = rho0 * g * dy * (Hwest**2 * ubarwest).sum()  # potential energy, sum in y
        PEwest[i] = rho0 * g * dy * (Hwest * ubarwest * zetawest).sum()  # potential energy anomaly, sum in y
        if kind == 'barotropic':
            KEwest[i] = 0.5 * rho0 * dy * (Hwest * ubarwest * (ubarwest**2 + vbarwest**2)).sum()  # sum in y
        elif kind == 'baroclinic':
            dzwest = Hwest/20.
            KEwest[i] = 0.5 * rho0 * dy * (dzwest * uwest * (uwest**2 + vwest**2)).sum()  # sum in z and y
        Pwest[i] = PEwest[i] + KEwest[i]

    for vec in teast:  # east
        ubareast = re(m['ubar'][0,vec[0]:vec[1],vec[2]:vec[2]+2],1).squeeze()  # interpolated to rho grid to match
        if kind == 'barotropic':
            vbareast = re(m['vbar'][0,vec[0]:vec[1]+1,vec[2]],0).squeeze()
        elif kind == 'baroclinic':
            ueast = re(m['u'][0,:, vec[0]:vec[1],vec[2]:vec[2]+2],2).squeeze()
            veast = re(m['v'][0,:, vec[0]:vec[1]+1,vec[2]],1).squeeze()
        zetaeast = m['zeta'][0,vec[0]:vec[1],vec[2]]
        Heast = (h[vec[0]:vec[1],vec[2]] + zetaeast)  # height
        # PEeast[i] = rho0 * g * dy * (Heast**2 * ubareast).sum()  # potential energy, sum in y
        PEeast[i] = rho0 * g * dy * (Heast * ubareast * zetaeast).sum()  # potential energy anomaly, sum in y
        if kind == 'barotropic':
            KEeast[i] = 0.5 * rho0 * dy * (Heast * ubareast * (ubareast**2 + vbareast**2)).sum()  # sum in y
        elif kind == 'baroclinic':
            dzeast = Heast/20.
            KEeast[i] = 0.5 * rho0 * dy * (dzeast * ueast * (ueast**2 + veast**2)).sum()  # sum in z and y
        Peast[i] = PEeast[i] + KEeast[i]

    for vec in tnorth:  # north
        vbarnorth = re(m['vbar'][0,vec[0]:vec[0]+2,vec[1]:vec[2]],0).squeeze()  # interpolated to rho grid to match
        if kind == 'barotropic':
            ubarnorth = re(m['ubar'][0,vec[0],vec[1]:vec[2]+1],0).squeeze()
        elif kind == 'baroclinic':
            unorth = re(m['u'][0,:, vec[0],vec[1]:vec[2]+1],1).squeeze()
            vnorth = re(m['v'][0,:, vec[0]:vec[0]+2,vec[1]:vec[2]],1).squeeze()
        zetanorth = m['zeta'][0,vec[0],vec[1]:vec[2]]
        Hnorth = (h[vec[0],vec[1]:vec[2]] + zetanorth)  # height
        # PEnorth[i] = rho0 * g * dy * (Hnorth**2 * vbarnorth).sum()  # potential energy, sum in y
        PEnorth[i] = rho0 * g * dx * (Hnorth * vbarnorth * zetanorth).sum()  # potential energy anomaly, sum in y
        if kind == 'barotropic':
            KEnorth[i] = 0.5 * rho0 * dx * (Hnorth * vbarnorth * (ubarnorth**2 + vbarnorth**2)).sum()  # sum in y
        elif kind == 'baroclinic':
            dznorth = Hnorth/20.
            KEnorth[i] = 0.5 * rho0 * dx * (dznorth * vnorth * (unorth**2 + vnorth**2)).sum()  # sum in z and y
        Pnorth[i] = PEnorth[i] + KEnorth[i]

    for vec in tsouth:  # south
        vbarsouth = re(m['vbar'][0,vec[0]:vec[0]+2,vec[1]:vec[2]],0).squeeze()  # interpolated to rho grid to match
        if kind == 'barotropic':
            ubarsouth = re(m['ubar'][0,vec[0],vec[1]:vec[2]+1],0).squeeze()
        elif kind == 'baroclinic':
            usouth = re(m['u'][0,:, vec[0],vec[1]:vec[2]+1],1).squeeze()
            vsouth = re(m['v'][0,:, vec[0]:vec[0]+2,vec[1]:vec[2]],1).squeeze()
        zetasouth = m['zeta'][0,vec[0],vec[1]:vec[2]]
        Hsouth = (h[vec[0],vec[1]:vec[2]] + zetasouth)  # height
        # PEsouth[i] = rho0 * g * dy * (Hsouth**2 * vbarsouth).sum()  # potential energy, sum in y
        PEsouth[i] = rho0 * g * dx * (Hsouth * vbarsouth * zetasouth).sum()  # potential energy anomaly, sum in y
        if kind == 'barotropic':
            KEsouth[i] = 0.5 * rho0 * dx * (Hsouth * vbarsouth * (ubarsouth**2 + vbarsouth**2)).sum()  # sum in y
        elif kind == 'baroclinic':
            dzsouth = Hsouth/20.
            KEsouth[i] = 0.5 * rho0 * dx * (dzsouth * vsouth * (usouth**2 + vsouth**2)).sum()  # sum in z and y
        Psouth[i] = PEsouth[i] + KEsouth[i]


    if doterms:

        ## Calculate mixing across inlet in time (buoyancy production)
        AKs = m['AKs'][0, 1:-1, 5:-4, 5:-4].squeeze()  # vertical viscosity for tracer
        rho = m['rho'][0, :, 5:-4, 5:-4].squeeze()  # density
        zeta = m['zeta'][0,5:-4,5:-4]
        dz = (h[5:-4,5:-4] + zeta)/20.  # easy since layers are uniform
        drhodz = (rho[:-1] - rho[1:])/dz
        Pmix[i] = g * dy * dx * (AKs*drhodz*dz).sum()
        ##

        ## Calculate momentum across inlet in time (shear production)
        AKv = m['AKv'][0, 1:-1, 5:-4, 5:-4].squeeze()  # vertical viscosity for momentum
        u = re(m['u'][0, :, 5:-4, 5:-2],2)  # to rho grid
        v = re(m['v'][0, :, 5:-2, 5:-4],1)  # to rho grid
        dudz = (u[:-1] - u[1:])/dz
        dvdz = (v[:-1] - v[1:])/dz
        Pmom[i] = rho0 * dy * dx * (AKv*(dudz**2 + dvdz**2)*dz).sum()
        ##

        ## Calculate power dissipation due to bottom friction, across domain
        # using velocity from bottom of water column
        Pfriction[i] = Cd * rho0 * dy * dx * abs((u[0]**2 + v[0]**2)**(3/2)).sum()
        ##

df = pd.DataFrame(index=dates)
df['PEeast'] = PEeast
df['KEeast'] = KEeast
df['Peast'] = Peast
df['PEwest'] = PEwest
df['KEwest'] = KEwest
df['Pwest'] = Pwest
df['PEnorth'] = PEnorth
df['KEnorth'] = KEnorth
df['Pnorth'] = Pnorth
df['PEsouth'] = PEsouth
df['KEsouth'] = KEsouth
df['Psouth'] = Psouth
df['Pin'] = Pwest - Pnorth
df['Pout'] = Psouth - Peast
df['P'] = df['Pin'] + df['Pout']
if doterms:
    df['Pmix'] = Pmix
    df['Pmom'] = Pmom
    df['Pfriction'] = Pfriction
df.to_csv('savedoutput/power_' + kind + which + '.csv')
# pd.read_csv('savedoutput/power.csv', index_col=0, parse_dates=True)

# np.savez('savedoutput/power_' + kind + '_temp.npz', Pwest=Pwest, Pnorth=Pnorth, Psouth=Psouth, dates=dates)
