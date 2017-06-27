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
import octant
import octant.depths
from plots import re


locgrid = '/pong/raid/kthyng/froude/ai65/grid.nc'
locmodel = '/pong/raid/kthyng/froude/ai65/OUT/'
grid = netCDF.Dataset(locgrid)
pm = grid['pm'][:]; pn = grid['pn'][:]
h = grid['h'][:]
lon_psi = grid['lon_psi'][:]; lat_psi = grid['lat_psi'][:]
lon_rho = grid['lon_rho'][:]; lat_rho = grid['lat_rho'][:]
dates = pd.date_range("2006-09-01", "2006-10-01", freq="15T")
ntimes = len(dates)
g = 9.81  # m/s^2
Cd = 3e-3  # friction coefficient
rho0 = 1023.7

Pin = np.empty(ntimes); Pout = np.empty(ntimes)
Pwest = np.empty((ntimes,20,146)); Pnorth = np.empty((ntimes,20,33))
Psouth = np.empty((ntimes,20,168))
Pmix = np.empty(ntimes); Pfriction = np.empty(ntimes)
Pmom = np.empty(ntimes)
for i in range(ntimes):

    m = netCDF.Dataset(locmodel + 'ocean_his_' + str(i+1).zfill(4) + '.nc')

    ## Calculate over power dissipation across inlet in time
    # west
    uwest = re(m['u'][0,:,296:-4,5:7],-1).squeeze()  # interpolated to rho grid to match
    zetawest = m['zeta'][0,296:-4,5]
    Wwest = (pn[296:-4,5]**-1)  # width
    Hwest = (h[296:-4,5] + zetawest)  # height
    # from Brian/Foreman
    Pwest[i] = (Wwest * Hwest * rho0 * uwest * (0.5*uwest**2 + g*zetawest))  # this is Lavelle's overall power
    # north
    vnorth = re(m['v'][0,:,-5:-3,5:38],-2).squeeze()
    zetanorth = m['zeta'][0,-5,5:38]
    Wnorth = (pm[-5,5:38]**-1)  # width
    Hnorth = (h[-5,5:38] + zetanorth)  # height
    Pnorth[i] = (Wnorth * Hnorth * rho0 * vnorth * (0.5*vnorth**2 + g*zetanorth))  # this is Lavelle's overall power
    # east is closed
    # south
    vsouth = re(m['v'][0,:,5:7,130:],-2).squeeze()
    zetasouth = m['zeta'][0,5,130:]
    Wsouth = (pm[5,130:]**-1)  # width
    Hsouth = (h[5,130:] + zetasouth)  # height
    Psouth[i] = (Wsouth * Hsouth * rho0 * vsouth * (0.5*vsouth**2 + g*zetasouth))  # this is Lavelle's overall power
    # Overall, then
    Pin[i] = Pwest[i].sum() - Pnorth[i].sum()  # west: + is in, for north: - is into domain
    Pout[i] = Psouth[i].sum()  # orientation makes this already opposite in sign
    ##

    ## Calculate mixing across inlet in time (buoyancy production)
    AKs = m['AKs'][0, 1:-1, 5:-4, 5:-4].squeeze()  # vertical viscosity for tracer
    rho = m['rho'][0, :, 5:-4, 5:-4].squeeze()  # density
    dz = h[5:-4,5:-4]/20.  # easy since layers are uniform
    drhodz = (rho[:-1] - rho[1:])*dz
    zeta = m['zeta'][0,5:-4,5:-4]
    A = (pm[5:-4,5:-4]**-1)*(pn[5:-4,5:-4]**-1)  # x-y area
    V = (zeta+h[5:-4,5:-4])*A  # volume, m^3
    Pmix[i] = (g*V*AKs*drhodz).sum()
    ##

    ## Calculate mixing across inlet in time (shear production)
    AKv = m['AKv'][0, 1:-1, 5:-4, 5:-4].squeeze()  # vertical viscosity for momentum
    u = re(m['u'][0, :, 5:-4, 5:-2],2)  # to rho grid
    v = re(m['v'][0, :, 5:-2, 5:-4],1)  # to rho grid
    dz = h[5:-4,5:-4]/20.  # easy since layers are uniform
    dudz = (u[:-1] - u[1:])*dz
    dvdz = (v[:-1] - v[1:])*dz
    zeta = m['zeta'][0,5:-4,5:-4]
    A = (pm[5:-4,5:-4]**-1)*(pn[5:-4,5:-4]**-1)  # x-y area
    V = (zeta + h[5:-4,5:-4])*A  # volume, m^3
    Pmom[i] = (rho0*V*AKv*(dudz**2 + dvdz**2)).sum()
    ##

    ## Calculate power dissipation due to bottom friction, across domain
    rhoxy = m['rho'][0, 0, 5:-4, 5:-4] + 1000  # density in bottom cell
    uxy = re(m['u'][0, 0, 5:-4, 5:-2],1)  # u in bottom cell, interpolated to rho grid
    vxy = re(m['v'][0, 0, 5:-2, 5:-4],0)  # v in bottom cell, interpolated to rho grid
    Pfriction[i] = (0.5*Cd*A*rhoxy*abs(np.sqrt(uxy**2 + vxy**2)**3)).sum()
    ##

df = pd.DataFrame(index=dates)
df['Pin'] = Pin
df['Pout'] = Pout
df['P'] = Pin + Pout
df['Pmix'] = Pmix
df['Pmom'] = Pmom
df['Pfriction'] = Pfriction
df.to_csv('savedoutput/power.csv')
# pd.read_csv('savedoutput/power.csv', index_col=0, parse_dates=True)

np.savez('savedoutput/power_temp.npz', Pwest=Pwest, Pnorth=Pnorth, Psouth=Psouth, dates=dates)
