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


# set up
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


def run(twest, teast, tnorth, tsouth, name, doterms=False):

    PEeast = np.zeros(ntimes); KEbteast = np.zeros(ntimes); KEbceast = np.zeros(ntimes)
    PEwest = np.zeros(ntimes); KEbtwest = np.zeros(ntimes); KEbcwest = np.zeros(ntimes)
    PEnorth = np.zeros(ntimes); KEbtnorth = np.zeros(ntimes); KEbcnorth = np.zeros(ntimes)
    PEsouth = np.zeros(ntimes); KEbtsouth = np.zeros(ntimes); KEbcsouth = np.zeros(ntimes)
    Pbtin = np.zeros(ntimes); Pbtout = np.zeros(ntimes)
    Pbcin = np.zeros(ntimes); Pbcout = np.zeros(ntimes)
    Pbtwest = np.zeros(ntimes); Pbtnorth = np.zeros(ntimes)
    Pbtsouth = np.zeros(ntimes); Pbteast = np.zeros(ntimes)
    Pbcwest = np.zeros(ntimes); Pbcnorth = np.zeros(ntimes)
    Pbcsouth = np.zeros(ntimes); Pbceast = np.zeros(ntimes)
    if doterms:
        Pmix = np.zeros(ntimes); Pfriction = np.zeros(ntimes)
        Pmom = np.zeros(ntimes)
    for i in range(ntimes):

        m = netCDF.Dataset(locmodel + 'ocean_his_' + str(i+1).zfill(4) + '.nc')
        # Calculate over power dissipation across inlet in time
        # Want to have all of these calculations on the psi grid to be at cell edges
        for vec in twest:  # west
            # in y direction need 2 extra indices: 1 for interpolating from rho to psi grid
            # and 1 for including the end point
            ubarwest = re(m['ubar'][0,vec[0]:vec[1]+2,vec[2]],0).squeeze()
            # in y direction need 1 extra index for including the end point
            # in x direction need 2 extra for interpolating from rho to psi grid
            vbarwest = re(m['vbar'][0,vec[0]:vec[1]+1,vec[2]:vec[2]+2],1).squeeze()
            uwest = re(m['u'][0,:, vec[0]:vec[1]+2,vec[2]],1).squeeze()
            vwest = re(m['v'][0,:, vec[0]:vec[1]+1,vec[2]:vec[2]+2],2).squeeze()
            # need to interpolate from rho to psi grid in both y and x
            zetawest = re(re(m['zeta'][0,vec[0]:vec[1]+2,vec[2]:vec[2]+2], 1), 0)
            Hwest = re(re(h[vec[0]:vec[1]+2,vec[2]:vec[2]+2], 1), 0) + zetawest  # height
            PEwest[i] = rho0 * g * dy * (Hwest * ubarwest * zetawest).sum()  # potential energy anomaly, sum in y
            KEbtwest[i] = 0.5 * rho0 * dy * (Hwest * ubarwest * (ubarwest**2 + vbarwest**2)).sum()  # sum in y
            dzwest = Hwest.squeeze()/20.
            KEbcwest[i] = 0.5 * rho0 * dy * (dzwest * uwest * (uwest**2 + vwest**2)).sum()  # sum in z and y
            Pbtwest[i] = PEwest[i] + KEbtwest[i]
            Pbcwest[i] = PEwest[i] + KEbcwest[i]

        for vec in teast:  # east
            ubareast = re(m['ubar'][0,vec[0]:vec[1]+2,vec[2]],0).squeeze()
            vbareast = re(m['vbar'][0,vec[0]:vec[1]+1,vec[2]:vec[2]+2],1).squeeze()
            ueast = re(m['u'][0,:, vec[0]:vec[1]+2,vec[2]],1).squeeze()
            veast = re(m['v'][0,:, vec[0]:vec[1]+1,vec[2]:vec[2]+2],2).squeeze()
            zetaeast = re(re(m['zeta'][0,vec[0]:vec[1]+2,vec[2]:vec[2]+2], 1), 0)
            Heast = re(re(h[vec[0]:vec[1]+2,vec[2]:vec[2]+2], 1), 0) + zetaeast  # height
            PEeast[i] = rho0 * g * dy * (Heast * ubareast * zetaeast).sum()  # potential energy anomaly, sum in y
            KEbteast[i] = 0.5 * rho0 * dy * (Heast * ubareast * (ubareast**2 + vbareast**2)).sum()  # sum in y
            dzeast = Heast.squeeze()/20.
            KEbceast[i] = 0.5 * rho0 * dy * (dzeast * ueast * (ueast**2 + veast**2)).sum()  # sum in z and y
            Pbteast[i] = PEeast[i] + KEbteast[i]
            Pbceast[i] = PEeast[i] + KEbceast[i]

        for vec in tnorth:  # north
            vbarnorth = re(m['vbar'][0,vec[0],vec[1]:vec[2]+2],0).squeeze()
            # in y direction 2 extra indices for intepolating from rho to psi grid with an extra index
            # in x direction 1 extra index to include end point
            ubarnorth = re(m['ubar'][0,vec[0]:vec[0]+2,vec[1]:vec[2]+1],0).squeeze()
            unorth = re(m['u'][0,:, vec[0]:vec[0]+2,vec[1]:vec[2]+1],1).squeeze()
            vnorth = re(m['v'][0,:, vec[0],vec[1]:vec[2]+2],1).squeeze()
            zetanorth = re(re(m['zeta'][0,vec[0]:vec[0]+2,vec[1]:vec[2]+2], 1), 0)
            Hnorth = re(re(h[vec[0]:vec[0]+2,vec[1]:vec[2]+2],1),0) + zetanorth  # height
            PEnorth[i] = rho0 * g * dx * (Hnorth * vbarnorth * zetanorth).sum()  # potential energy anomaly, sum in y
            KEbtnorth[i] = 0.5 * rho0 * dx * (Hnorth * vbarnorth * (ubarnorth**2 + vbarnorth**2)).sum()  # sum in y
            dznorth = Hnorth.squeeze()/20.
            KEbcnorth[i] = 0.5 * rho0 * dx * (dznorth * vnorth * (unorth**2 + vnorth**2)).sum()  # sum in z and y
            Pbtnorth[i] = PEnorth[i] + KEbtnorth[i]
            Pbcnorth[i] = PEnorth[i] + KEbcnorth[i]

        for vec in tsouth:  # south
            vbarsouth = re(m['vbar'][0,vec[0],vec[1]:vec[2]+2],0).squeeze()
            ubarsouth = re(m['ubar'][0,vec[0]:vec[0]+2,vec[1]:vec[2]+1],0).squeeze()
            usouth = re(m['u'][0,:, vec[0]:vec[0]+2,vec[1]:vec[2]+1],1).squeeze()
            vsouth = re(m['v'][0,:, vec[0],vec[1]:vec[2]+2],1).squeeze()
            zetasouth = re(re(m['zeta'][0,vec[0]:vec[0]+2,vec[1]:vec[2]+2], 1), 0)
            Hsouth = re(re(h[vec[0]:vec[0]+2,vec[1]:vec[2]+2],1),0) + zetasouth  # height
            PEsouth[i] = rho0 * g * dx * (Hsouth * vbarsouth * zetasouth).sum()  # potential energy anomaly, sum in y
            KEbtsouth[i] = 0.5 * rho0 * dx * (Hsouth * vbarsouth * (ubarsouth**2 + vbarsouth**2)).sum()  # sum in y
            dzsouth = Hsouth.squeeze()/20.
            KEbcsouth[i] = 0.5 * rho0 * dx * (dzsouth * vsouth * (usouth**2 + vsouth**2)).sum()  # sum in z and y
            Pbtsouth[i] = PEsouth[i] + KEbtsouth[i]
            Pbcsouth[i] = PEsouth[i] + KEbcsouth[i]

        if doterms:

            ## Calculate mixing across inlet in time (buoyancy production)
            # in y direction +1 index at south to include rho inside of psi edge
            # and +1 index at north to include end point. (Same for x)
            AKs = m['AKs'][0, 1:-1, tsouth[0][0]+1:tnorth[0][0]+1, twest[0][2]+1:teast[0][2]+1]  # vertical viscosity for tracer
            rho = m['rho'][0, :, tsouth[0][0]+1:tnorth[0][0]+1, twest[0][2]+1:teast[0][2]+1]  # density
            zeta = m['zeta'][0,tsouth[0][0]+1:tnorth[0][0]+1, twest[0][2]+1:teast[0][2]+1]
            dz = (h[tsouth[0][0]+1:tnorth[0][0]+1, twest[0][2]+1:teast[0][2]+1] + zeta)/20.  # easy since layers are uniform
            drhodz = (rho[:-1] - rho[1:])/dz
            Pmix[i] = g * dy * dx * (AKs*drhodz*dz).sum()
            ##

            ## Calculate momentum across inlet in time (shear production)
            AKv = m['AKv'][0, 1:-1, tsouth[0][0]+1:tnorth[0][0]+1, twest[0][2]+1:teast[0][2]+1]  # vertical viscosity for momentum
            u = re(m['u'][0, :, tsouth[0][0]+1:tnorth[0][0]+1, twest[0][2]:teast[0][2]+1],2)
            v = re(m['v'][0, :, tsouth[0][0]:tnorth[0][0]+1, twest[0][2]+1:teast[0][2]+1],1)
            dudz = (u[:-1] - u[1:])/dz
            dvdz = (v[:-1] - v[1:])/dz
            Pmom[i] = rho0 * dy * dx * (AKv*(dudz**2 + dvdz**2)*dz).sum()
            ##

            ## Calculate power dissipation due to bottom friction, across domain
            # using velocity from bottom of water column
            Pfriction[i] = Cd * rho0 * dy * dx * abs((u[0]**2 + v[0]**2)**(3/2)).sum()
            ##

    df = pd.DataFrame(index=dates)
    df['PEeast'] = PEeast; df['KEbteast'] = KEbteast; df['KEbceast'] = KEbceast
    df['Pbteast'] = Pbteast; df['Pbceast'] = Pbceast;
    df['PEwest'] = PEwest; df['KEbtwest'] = KEbtwest; df['KEbcwest'] = KEbcwest
    df['Pbtwest'] = Pbtwest; df['Pbcwest'] = Pbcwest;
    df['PEnorth'] = PEnorth; df['KEbtnorth'] = KEbtnorth; df['KEbcnorth'] = KEbcnorth
    df['Pbtnorth'] = Pbtnorth; df['Pbcnorth'] = Pbcnorth;
    df['PEsouth'] = PEsouth; df['KEbtsouth'] = KEbtsouth; df['KEbcsouth'] = KEbcsouth
    df['Pbtsouth'] = Pbtsouth; df['Pbcsouth'] = Pbcsouth;
    if np.isnan(Pbtwest.sum()):
        Pbtwest = np.zeros(Pbtwest.size)
        Pbcwest = np.zeros(Pbcwest.size)
    if np.isnan(Pbteast.sum()):
        Pbteast = np.zeros(Pbteast.size)
        Pbceast = np.zeros(Pbceast.size)
    if np.isnan(Pbtnorth.sum()):
        Pbtnorth = np.zeros(Pbtnorth.size)
        Pbcnorth = np.zeros(Pbcnorth.size)
    if np.isnan(Pbtsouth.sum()):
        Pbtsouth = np.zeros(Pbtsouth.size)
        Pbcsouth = np.zeros(Pbcsouth.size)
    df['Pbtin'] = Pbtwest - Pbtnorth
    df['Pbcin'] = Pbcwest - Pbcnorth
    df['Pbtout'] = Pbtsouth - Pbteast
    df['Pbcout'] = Pbcsouth - Pbceast
    df['Pbt'] = df['Pbtin'] + df['Pbtout']
    df['Pbc'] = df['Pbcin'] + df['Pbcout']
    if doterms:
        df['Pmix'] = Pmix
        df['Pmom'] = Pmom
        df['Pfriction'] = Pfriction
    df.to_csv('savedoutput/power/' + name + '.csv')
    # pd.read_csv('savedoutput/power.csv', index_col=0, parse_dates=True)


if __name__ == "__main__":


    # parse input arguments
    which = 'stripsd50' # 'all' 'center'
    doterms = True


    if which == 'testtop':
        # west transects: [0:1, 2]
        twest = np.array([[335, 375, 5]])
        # east transects: [0:1, 2]
        teast = np.array([[335, 375, 280]])
        # north transects: [0, 1:2]
        tnorth = np.array([[375,5,280]])
        # south transects: [0, 1:2]
        tsouth = np.array([[335, 5,280]])
    elif which == 'testbottom':
        # west transects: [0:1, 2]
        twest = np.array([[295, 335, 5]])
        # east transects: [0:1, 2]
        teast = np.array([[295, 335, 280]])
        # north transects: [0, 1:2]
        tnorth = np.array([[335,5,280]])
        # south transects: [0, 1:2]
        tsouth = np.array([[295, 5,280]])
    elif which == 'all':
        twest = np.array([[0, 444, 0]])
        teast = np.array([[0, 444, 296]])
        tnorth = np.array([[444, 0, 296]])
        tsouth = np.array([[0, 0, 296]])



    # choose domain to use, possibly looping. Run function.
    if 'strips' not in which:
        run(twest, teast, tnorth, tsouth, which)

    elif which == 'strips':

        for j in range(lon_psi.shape[0]):
            name = which + str(j)
            fname = 'savedoutput/power/' + name + '.csv'
            if os.path.exists(fname):
                continue
            twest = np.array([[j, j+1, 0]])  # west transects: [0:1, 2]
            teast = np.array([[j, j+1, 296]])  # east transects: [0:1, 2]
            tnorth = np.array([[j+1, 0, 296]])  # north transects: [0, 1:2]
            tsouth = np.array([[j, 0, 296]])  # south transects: [0, 1:2]
            run(twest, teast, tnorth, tsouth, name, doterms=doterms)

    elif which == 'stripsd10':
        dd = 10
        for j in range(lon_psi.shape[0])[::dd]:
            twest = np.array([[j, j+dd, 0]])  # west transects: [0:1, 2]
            teast = np.array([[j, j+dd, 296]])  # east transects: [0:1, 2]
            tnorth = np.array([[j+dd, 0, 296]])  # north transects: [0, 1:2]
            tsouth = np.array([[j, 0, 296]])  # south transects: [0, 1:2]
            run(twest, teast, tnorth, tsouth, which + str(j), doterms=doterms)

    elif which == 'stripsd20':
        dd = 20
        for j in range(lon_psi.shape[0])[::dd]:
            twest = np.array([[j, j+dd, 0]])  # west transects: [0:1, 2]
            teast = np.array([[j, j+dd, 296]])  # east transects: [0:1, 2]
            tnorth = np.array([[j+dd, 0, 296]])  # north transects: [0, 1:2]
            tsouth = np.array([[j, 0, 296]])  # south transects: [0, 1:2]
            run(twest, teast, tnorth, tsouth, which + str(j), doterms=doterms)

    elif which == 'stripsd50':
        dd = 50
        for j in range(lon_psi.shape[0])[::dd]:
            twest = np.array([[j, j+dd, 0]])  # west transects: [0:1, 2]
            teast = np.array([[j, j+dd, 296]])  # east transects: [0:1, 2]
            tnorth = np.array([[j+dd, 0, 296]])  # north transects: [0, 1:2]
            tsouth = np.array([[j, 0, 296]])  # south transects: [0, 1:2]
            run(twest, teast, tnorth, tsouth, which + str(j), doterms=doterms)
