'''
Plotting.
Example:
python3 plots.py 'vort' 'full' 'bar' 'intime'
python3 plots.py 'mix' 'full' 'bar' 'intime'
python3 plots.py 'mom' 'full' 'bar' 'intime'
python3 plots.py 'friction' 'full' 'bar' 'intime'
python3 plots.py 'rho' 'full' 'surface' 'intime'
python3 plots.py 'upwelling' 'full' 'bar' 'intime'
python3 plots.py 'upsloping' 'full' 'bar' 'intime'
python3 plots.py 'mix' 'full' 'bar' 'mean'
python3 plots.py 'upwelling' 'full' 'bar' 'mean'
python3 plots.py 'upsloping' 'full' 'bar' 'mean'
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

g = 9.81  # m/s^2
Cd = 3e-3  # friction coefficient
rho0 = 1023.7

merc = ccrs.Mercator()
pc = ccrs.PlateCarree()
locgrid = '/pong/raid/kthyng/froude/ai65/grid.nc'
locmodel = '/pong/raid/kthyng/froude/ai65/OUT/'
try:
    grid = netCDF.Dataset(locgrid)
except:
    grid = netCDF.Dataset('grid.nc')
pm = grid['pm'][0,0]; pn = grid['pn'][0,0]
dx = pm**-1; dy = pn**-1

lon_psi = grid['lon_psi'][:]; lat_psi = grid['lat_psi'][:]
lon_rho = grid['lon_rho'][:]; lat_rho = grid['lat_rho'][:]
dates = pd.date_range("2006-09-01", "2006-10-01", freq="15T")
h = grid['h'][:]
ind = (lon_rho > -122.64)*(lat_rho > 48.17)
h[ind] = np.nan
ind = (lon_rho > -122.58)*(lat_rho > 48.0)
h[ind] = np.nan
def get_zeta():
    '''Read in zeta time series for plotting.'''
    # lon_zeta = -122.71497166666667; lat_zeta = 48.149936666666669
    # middle of channel

    fname = 'savedoutput/tidal.csv'
    if os.path.exists(fname):
        df = pd.read_csv(fname, index_col=0, parse_dates=True)
    else:
        # read in tidal signal for plotting
        j, i = 313, 98
        zeta = np.empty(len(dates))
        for ind in range(len(dates)):
            mzeta = netCDF.Dataset(locmodel + 'ocean_his_' + str(ind+1).zfill(4) + '.nc')
            zeta[ind] = mzeta['zeta'][0, j, i]
            mzeta.close()
        df = pd.DataFrame(index=dates)
        df['zeta'] = zeta
        df.to_csv(fname)
    return df

def re(A, dim):
    """
    Average neighboring elements in an array A along a dimension dim.
    Args:
        A (array): array size, [m x n] in 2D case. Can be up to 3D.
        dim (int): dimension on which to act. Can be up to 2 (0, 1, 2).
    Returns:
        * A - array of size [(m-1) x n] if dim is 0
    """

    # B is A but rolled such that the dimension that is to be resized is in
    # the 0 position
    B = np.rollaxis(A, dim)
    # Do averaging
    B = 0.5*(B[0:-1]+B[1:])
    # Roll back to original
    return np.rollaxis(B, 0, dim+1)


def setup(zoom='full', dobathy=False):
    '''Setup plot.'''

    if zoom == 'full':
        figsize = 5.25, 7
        box = [-122.8, -122.54, 47.9655, 48.227]
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(top=0.96, bottom=0.02, left=0.06, right=0.99)
    ax = fig.add_subplot(111, projection=merc)#, facecolor='0.8')
    ax.set_extent(box, pc)
    gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
    gl.xlocator = mticker.FixedLocator([-122.8, -122.75, -122.7, -122.65, -122.6, -122.55, -122.5])
    # the following two make the labels look like lat/lon format
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabels_top = False  # turn off labels where you don't want them
    gl.ylabels_right = False
    # plot coastline
    ax.contour(grid['lon_rho'][:], grid['lat_rho'][:], grid['mask_rho'][:], [0], colors='0.3', transform=pc, linewidths=1.5)
    # plot bathymetry
    if dobathy:
        ax.contour(grid['lon_rho'][:], grid['lat_rho'][:], h, np.arange(25,200,25), colors='0.3', transform=pc, linewidths=0.25)

    return fig, ax


def setupvar(varname, kind):
    '''Setup variable and other choices given the variable plotted.'''

    if varname == 'vort':
        cmap = cmo.curl
        vmax = 0.01
        vmin = -vmax
        lon = lon_rho
        lat = lat_rho
        label = 'Vertical\nvorticity\n[s$^{-1}$]'
        ticks = np.linspace(vmin, vmax, 5)
        ctidal = '#86174F' ## pink or '#107469' teal
        alpha = 1
    elif varname == 'mix':
        cmap = cmo.amp
        if kind == 'intime':
            vmax = 2.0
        elif kind == 'mean':
            vmax = 0.5
        vmin = 0
        lon = lon_psi
        lat = lat_psi
        label = 'Dissipation:\nbuoyancy\nproduction\n[W/m$^2$]'
        # label = 'Vertical\nmixing\n[kg$\cdot$m$^{-2}\cdot$s$^{-1}$]'
        ticks = np.linspace(vmin, vmax, 6)
        ctidal = '#71001D'
        alpha = 0.6
    elif varname == 'mom':
        cmap = cmo.amp
        if kind == 'intime':
            vmax = 5
        elif kind == 'mean':
            vmax = 1.25
        vmin = 0
        lon = lon_psi
        lat = lat_psi
        label = 'Dissipation:\nshear\nproduction\n[W/m$^2$]'
        # label = 'Vertical\nmixing\n[kg$\cdot$m$^{-2}\cdot$s$^{-1}$]'
        ticks = np.linspace(vmin, vmax, 6)
        ctidal = '#71001D'
        alpha = 0.6
    elif varname == 'friction':
        cmap = cmo.amp
        if kind == 'intime':
            vmax = 8.0
        elif kind == 'mean':
            vmax = 2.0
        vmin = 0
        lon = lon_psi
        lat = lat_psi
        label = 'Dissipation:\nbottom\nfriction\n[W/m$^2$]'
        # label = 'Vertical\nmixing\n[kg$\cdot$m$^{-2}\cdot$s$^{-1}$]'
        ticks = np.linspace(vmin, vmax, 6)
        ctidal = '#71001D'
        alpha = 0.6
    elif varname == 'rho':
        cmap = cmo.haline
        vmax = 24.5
        vmin = 21
        lon = lon_psi
        lat = lat_psi
        label = 'Density\n[kg$\cdot$ m$^{-3}$]'
        ticks = np.linspace(vmin, vmax, 6)
        ctidal = '#0A3885'
        alpha = 0.6
    elif varname == 'upwelling':
        cmap = cmo.delta
        if kind == 'intime':
            vmax = 0.1
        elif kind == 'mean':
            vmax = 0.05
        vmin = -vmax
        lon = lon_psi
        lat = lat_psi
        label = 'Upwelling\nvelocity\n[m$\cdot$s$^{-1}$]'
        ticks = np.linspace(vmin, vmax, 6)
        ctidal = '#226E13'
        alpha = 0.6
    elif varname == 'upsloping':
        cmap = cmo.delta
        if kind == 'intime':
            vmax = 0.1
        elif kind == 'mean':
            vmax = 0.05
        vmin = -vmax
        lon = lon_psi
        lat = lat_psi
        label = 'Upsloping\nvelocity\n[m$\cdot$s$^{-1}$]'
        ticks = np.linspace(vmin, vmax, 6)
        ctidal = '#226E13'
        alpha = 0.6

    return lon, lat, cmap, vmin, vmax, label, ticks, ctidal, alpha


def setupvarforstep(m, varname, depth='bar'):
    '''Setup variable and other choices given the variable plotted.'''

    if varname == 'vort':
        if depth == 'bar':
            v = m['vbar'][:].squeeze()
            u = m['ubar'][:].squeeze()
        elif depth == 'surface':
            v = m['v'][-1,:,:]
        var = (v[:,:-1] - v[:,1:])*re(re(pm,1),0) - (u[:-1,:] - u[1:,:])*re(re(pn,1),0)

    elif varname == 'mix':
        if depth == 'bar':
            AKs = m['AKs'][0, 1:-1, 1:-1, 1:-1].squeeze()  # vertical viscosity for scalar
            rho = m['rho'][0, :, 1:-1, 1:-1].squeeze()  # salt
            zeta = m['zeta'][0, 1:-1, 1:-1].squeeze()
            dz = (h[1:-1,1:-1] + zeta)/20.  # easy since layers are uniform
            dsdz = (rho[:-1] - rho[1:])/dz
            var = g * (AKs * dsdz * dz).sum(axis=0)  # sum in z

    elif varname == 'mom':
        if depth == 'bar':
            AKv = m['AKv'][0, 1:-1, 1:-1, 1:-1].squeeze()  # vertical viscosity for momentum
            u = re(m['u'][0, :, 1:-1, :],2).squeeze()
            v = re(m['v'][0, :, :, 1:-1],1).squeeze()
            zeta = m['zeta'][0, 1:-1, 1:-1].squeeze()
            dz = (h[1:-1,1:-1] + zeta)/20.  # easy since layers are uniform
            dudz = (u[:-1] - u[1:])/dz
            dvdz = (v[:-1] - v[1:])/dz
            var = rho0 * (AKv * (dudz**2 + dvdz**2) * dz).sum(axis=0)  # sum in z

    elif varname == 'friction':
        if depth == 'bar':  # not really depth averaged
            u = re(m['u'][0, 0, 1:-1, :],1).squeeze()
            v = re(m['v'][0, 0, :, 1:-1],0).squeeze()
            var = Cd * rho0 * abs((u**2 + v**2)**(3/2))

    elif varname == 'rho':
        # if depth == 'bar':
        #     v = m['vbar'][:].squeeze()
        #     u = m['ubar'][:].squeeze()
        if depth == 'surface':
            # Note: skip ghost cells in x and y so that can properly plot grid cell boxes with pcolormesh
            var = m['rho'][0,-1,1:-1,1:-1].squeeze()

    elif varname == 'upwelling':
        if depth == 'bar':
            var = m['omega'][0,:,1:-1,1:-1].squeeze().mean(axis=0)
        # if depth == 'surface':
        #     # Note: skip ghost cells in x and y so that can properly plot grid cell boxes with pcolormesh
        #     var = m['rho'][0,-1,1:-1,1:-1].squeeze()

    elif varname == 'upsloping':
        if depth == 'bar':
            w = m['w'][0,:,1:-1,1:-1].squeeze()
            uw = m['omega'][0,:,1:-1,1:-1].squeeze()
            var = (w - uw).mean(axis=0)
        # if depth == 'surface':
        #     # Note: skip ghost cells in x and y so that can properly plot grid cell boxes with pcolormesh
        #     var = m['rho'][0,-1,1:-1,1:-1].squeeze()

    return var


def add_stuff(lon, lat, vart, cmap, vmin, vmax, label, ticks, alpha):
    '''Plot actual stuff.'''

    mapp = ax.pcolormesh(lon, lat, vart, cmap=cmap, vmin=vmin, vmax=vmax, transform=pc)

    # Colorbar in lower left corner
    cax = fig.add_axes([0.795, 0.25, 0.02, 0.35]) #colorbar axes
    # cax = fig.add_axes([0.09, 0.15, 0.25, 0.015]) #colorbar axes
    cb = fig.colorbar(mapp, cax=cax, orientation='vertical')
    # cb.set_label(label, fontsize=12, color='0.2')
    cb.ax.tick_params(labelsize=10, length=2, color='0.2', labelcolor='0.2')
    cb.set_ticks(ticks)
    # change colorbar tick color http://stackoverflow.com/questions/9662995/matplotlib-change-title-and-colorbar-text-and-tick-colors
    cbtick = plt.getp(cb.ax.axes, 'yticklabels')
    plt.setp(cbtick, color='0.2')
    # text above colorbar
    cax.text(1.5, 1.1, label, ha='center', fontsize=12, transform=cax.transAxes)


def add_tidal(date, ctidal):
    '''Add tidal signal.'''

    axt = fig.add_axes([0.43, 0.85, 0.5, .11], facecolor='None')
    # turn frame off
    for axis in ['top','left','right','bottom']:
        axt.spines[axis].set_linewidth(0.0)
    # control y ticks
    axt.set_yticks([-1.5, 0])  # max is skewed off top of plot
    axt.tick_params(axis='y', pad=0, labelsize=8, length=0, color='0.5', labelcolor='0.5')
    # plot tide signal
    df.plot(ax=axt, color='0.5', lw=1, legend=False, rot=-30)
    df[date:date].plot(marker='o', ax=axt, ms=4, color=ctidal, legend=False)
    axt.xaxis.set_major_locator(Dates.DayLocator())
    axt.xaxis.set_major_formatter(Dates.DateFormatter('%d'))
    axt.text(0.55, 0.8, 'September 2006', transform=ax.transAxes, color='0')
    fig.text(0.12, 0.025, date, transform=ax.transAxes, color='0', fontsize=12)


def add_arrows(ax, m, depth='bar', zoom='full'):
    '''add current quiver arrows.'''

    if zoom == 'full':
        ddx = 8; ddy = 8
        qkx, qky = 0.15, 0.08  # quiver key location

    if depth == 'bar':
        v = m['vbar'][:].squeeze()
        u = m['ubar'][:].squeeze()
    elif depth == 'surface':
        v = m['v'][0,-1,:,:].squeeze()
        u = m['u'][0,-1,:,:].squeeze()
    u, v = re(u,0), re(v,1)  # resize to psi grid

    Q = ax.quiver(lon_psi[::ddy, ::ddx], lat_psi[::ddy, ::ddx],
                  u[::ddy, ::ddx], v[::ddy, ::ddx], color='k', alpha=0.3,
                  pivot='middle', scale=25, transform=pc)
    qk = ax.quiverkey(Q, qkx, qky, 1.0, r'1.0 m$\cdot$s$^{-1}$ current',
                      labelcolor='0.5', fontproperties={'size': '10'})


if __name__ == "__main__":

    # parse the input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('var', type=str, help='variable to plot', default="salt")
    parser.add_argument('zoom', type=str, help='which zoom to use: "full"', default="full")
    parser.add_argument('depth', type=str, help='which depth: "bar", "surface"', default="bar")
    parser.add_argument('kind', type=str, help='which kind over time: "intime", "mean" across 30 days', default="intime")
    # parser.add_argument('--dend', type=str, help='dend', default=None)
    args = parser.parse_args()
    varname = args.var
    zoom = args.zoom
    depth = args.depth
    kind = args.kind

    # prep variable
    lon, lat, cmap, vmin, vmax, label, ticks, ctidal, alpha = setupvar(varname, kind)
    # tidal signal
    df = get_zeta()

    # make directories
    basename = 'figures/' + varname + '/' + zoom + '/' + depth + '/'
    basenames = basename.split('/')
    for i in range(len(basenames)-1):
        tempname = '/'.join(basenames[:i+1])
        if not os.path.exists(tempname):
            os.makedirs(tempname)

    if varname == 'upwelling' or varname == 'upsloping':
        dobathy = True
    else:
        dobathy = False

    ntimes = len(dates)

    if kind == 'mean':
        Var = np.zeros((lon.shape[0]-1, lon.shape[1]-1))
    for i in range(ntimes):
        if kind == 'intime':
            date = pd.to_datetime(dates[i]).strftime('%b %02d, %Y %H:%M')
            fname = basename + pd.to_datetime(dates[i]).strftime('%Y-%m-%02dT%H:%M')
        elif kind == 'mean':
            fname = basename + 'mean'
        if os.path.exists(fname + '.png'):
            continue

        m = netCDF.Dataset(locmodel + 'ocean_his_' + str(i+1).zfill(4) + '.nc')
        vart = setupvarforstep(m, varname, depth=depth)
        if kind == 'mean':
            meanname = 'savedoutput/' + varname + '_mean.npz'
            if os.path.exists(meanname):
                i = ntimes-1  # short-circuit loop and jump to end
                meanfile = np.load(meanname)
                vart = meanfile['var']
            else:
                Var += vart
                # on the last time step, complete average and rename
                if i == ntimes-1:
                    vart = Var/ntimes
                    np.savez(meanname, lon=lon, lat=lat, var=vart)
                else:
                    continue

        fig, ax = setup(zoom=zoom, dobathy=dobathy)
        add_stuff(lon, lat, vart, cmap, vmin, vmax, label, ticks, alpha)
        if kind == 'intime':
            add_tidal(date, ctidal) # need to save a tidal signal to use
            add_arrows(ax, m, depth=depth, zoom=zoom)
        fig.savefig(fname + '.png', bbox_inches='tight', dpi=100)
        plt.close(fig)
