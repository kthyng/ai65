'''
Plotting.
Example:
python plots.py 'vort' 'full'
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


merc = ccrs.Mercator()
pc = ccrs.PlateCarree()
locgrid = '/pong/raid/kthyng/froude/ai65/grid.nc'
locmodel = '/pong/raid/kthyng/froude/ai65/OUT/'
grid = netCDF.Dataset(locgrid)
pm = grid['pm'][:]; pn = grid['pn'][:]
lon_psi = grid['lon_psi'][:]; lat_psi = grid['lat_psi'][:]
lon_rho = grid['lon_rho'][:]; lat_rho = grid['lat_rho'][:]
land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])

# model tidal signal for later
# read in saved zeta
# mat = scipy.io.loadmat('savedoutput/zeta_adcpzetasites.mat')
# location 0 is in the middle of the channel between Admiralty Head and Point Wilson
dates = pd.date_range("2006-09-01", "2006-10-01", freq="15T")
# df = pd.DataFrame(index=dates)
# df['zeta'] = mat['data'][:,0]
# zeta = mat['data'][:,0]
# zmax = max(abs(zeta))

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


def setup(zoom='full'):
    '''Setup plot.'''

    if zoom == 'full':
        figsize = 5.25, 7
        box = [-122.8, -122.54, 47.9655, 48.227]
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(top=0.96, bottom=0.04, left=0.08, right=0.99)
    ax = fig.add_subplot(111, projection=merc)#, facecolor='0.8')
    # ax.set_facecolor("0.8")
    # ax.background_patch.set_fill(True)
    # Set RGB value to land color
    # ax.imshow(np.tile(np.array([[[235, 235, 235]]],
    #           dtype=np.uint8), [2, 2, 1]),
    #       origin='center',
    #       transform=pc,
    #       extent=box)

    ax.set_extent(box, pc)
    gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
    gl.xlocator = mticker.FixedLocator([-122.8, -122.75, -122.7, -122.65, -122.6, -122.55, -122.5])
    # the following two make the labels look like lat/lon format
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabels_top = False  # turn off labels where you don't want them
    gl.ylabels_right = False

    # plot isobaths
    # ax.contour(grid['lon_rho'][:], grid['lat_rho'][:], grid['h'][:], np.arange(0, 200, 20), colors='0.6', transform=pc, linewidths=0.5)

    # plot coastline
    ax.contour(grid['lon_rho'][:], grid['lat_rho'][:], grid['mask_rho'][:], [0], colors='0.3', transform=pc, linewidths=2, alpha=0.7)

    return fig, ax


def setupvar(varname):
    '''Setup variable and other choices given the variable plotted.'''

    if varname == 'vort':
        cmap = cmo.curl
        vmax = 0.01
        vmin = -vmax
        lon = lon_rho
        lat = lat_rho
        label = r'Vertical vorticity [s$^{-1}$]'
        ticks = np.linspace(vmin, vmax, 5)
        ctidal = '#86174F' ## pink or '#107469' teal

    # elif varname == 'salt':

    return lon, lat, cmap, vmin, vmax, label, ticks, ctidal


def setupvarforstep(m, varname, depth='bar'):
    '''Setup variable and other choices given the variable plotted.'''

    if varname == 'vort':
        if depth == 'bar':
            v = m['vbar'][:].squeeze()
            u = m['ubar'][:].squeeze()
        elif depth == 'surface':
            v = m['v'][-1,:,:]
        var = (v[:,:-1] - v[:,1:])*re(re(pm,1),0) - (u[:-1,:] - u[1:,:])*re(re(pn,1),0)

    # elif varname == 'mix'
    # elif varname == 'salt':
        # if depth == 'bar':
        #     v = m['v'].sel(ocean_time=date)
        # elif depth == 'surface':
        #     v = m['v'].sel(ocean_time=date).isel(s_rho=-1, eta_rho=slice(1,-1), xi_rho=slice(1,-1))
    # Note: skip ghost cells in x and y so that can properly plot grid cell boxes with pcolormesh
    # salt = m.salt.sel(ocean_time=plotdate).isel(s_rho=-1, eta_rho=slice(1,-1), xi_rho=slice(1,-1))

    return var


def add_stuff(lon, lat, vart, cmap, vmin, vmax, label, ticks):
    '''Plot actual stuff.'''

    mapp = ax.pcolormesh(lon, lat, vart, cmap=cmap, vmin=vmin, vmax=vmax, transform=pc)
    # ax.add_feature(land_10m, facecolor='0.8')
    # ax.coastlines(resolution='10m')  # coastline resolution options are '110m', '50m', '10m'
    # land data from: http://openstreetmapdata.com/data/land-polygons
    # example from: https://ocefpaf.github.io/python4oceanographers/blog/2015/06/22/osm/
    # shp = shapereader.Reader('data/land-polygons-complete-4326/land_polygons')
    # for record, geometry in zip(shp.records(), shp.geometries()):
    #     ax.add_geometries([geometry], pc, facecolor='0.8', edgecolor='black')
    # # coastline
    # shp = shapereader.Reader('./data/OSM_coastline/BTS')
    # for record, geometry in zip(shp.records(), shp.geometries()):
    #     ax.add_geometries([geometry], ccrs.PlateCarree(), facecolor='w',
    #                       edgecolor='black')

    ax.set_title(label)

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


def add_tidal(date, ctidal):
    '''Add tidal signal.'''

    axt = fig.add_axes([0.43, 0.85, 0.5, .11], facecolor='None')
    # axt.background_patch.set_fill(False)
    # turn frame off
    for axis in ['top','left','right','bottom']:
        axt.spines[axis].set_linewidth(0.0)
        # axr.spines[axis].set_linewidth(0.05)
    # axr.spines['bottom'].set_linewidth(0.0)
    # ticks

    # axt.yaxis.tick_right()
    axt.set_yticks([-1.5, 0])  # max is skewed off top of plot
    axt.tick_params(axis='y', pad=0, labelsize=8, length=0, color='0.5', labelcolor='0.5')

    # axt.get_yaxis().set_visible(False)
    # axt.get_xaxis().set_visible(False)
    df.plot(ax=axt, color='0.5', lw=0.9, legend=False, rot=-30)
    axt.plot(date, df[date:date], 'o', ms=5, color=ctidal)
    axt.xaxis.set_major_locator(Dates.DayLocator())
    axt.xaxis.set_major_formatter(Dates.DateFormatter('%d'))
    axt.text(0.55, 0.8, 'September 2006', transform=ax.transAxes)
    # # the plot itself, in indices
    # axt.plot(np.arange(zeta.size), zeta, '0.5', lw=0.9)
    # axt.plot(itime, zeta[itime], 'o', ms=5, color=ctidal)


def add_arrows(ax, m, depth='bar', zoom='full'):
    '''add current quiver arrows.'''

    if zoom == 'full':
        ddx = 12; ddy = 12
        qkx, qky = 0.17, 0.04  # quiver key location

    if depth == 'bar':
        v = m['vbar'][:].squeeze()
        u = m['ubar'][:].squeeze()
    elif depth == 'surface':
        v = m['v'][-1,:,:]
        u = m['u'][-1,:,:]
    u, v = re(u,0), re(v,1)  # resize to psi grid

    Q = ax.quiver(lon_psi[::ddy, ::ddx], lat_psi[::ddy, ::ddx],
                  u[::ddy, ::ddx], v[::ddy, ::ddx], color='k', alpha=0.4,
                  pivot='middle', scale=20, transform=pc)
    qk = ax.quiverkey(Q, qkx, qky, 0.5, r'0.5 m$\cdot$s$^{-1}$ current',
                      labelcolor='0.5', fontproperties={'size': '10'})




if __name__ == "__main__":

    # parse the input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('var', type=str, help='variable to plot', default="salt")
    parser.add_argument('zoom', type=str, help='which zoom to use: "full"', default="full")
    parser.add_argument('depth', type=str, help='which depth: "bar", "surface"', default="bar")
    # parser.add_argument('--dend', type=str, help='dend', default=None)
    args = parser.parse_args()
    varname = args.var
    zoom = args.zoom
    depth = args.depth

    # prep variable
    lon, lat, cmap, vmin, vmax, label, ticks, ctidal = setupvar(varname)

    df = get_zeta()

    fig, ax = setup(zoom=zoom)
    # for i, date in zip(range(2881), dates):
    i=1000; date = pd.to_datetime(dates[i]).strftime('%Y-%m-%02d %H:%M')
    m = netCDF.Dataset(locmodel + 'ocean_his_' + str(i+1).zfill(4) + '.nc')
    # dates = m['ocean_time'].data
    vart = setupvarforstep(m, varname, depth=depth)
    add_stuff(lon, lat, vart, cmap, vmin, vmax, label, ticks)
    add_tidal(date, ctidal) # need to save a tidal signal to use
    add_arrows(ax, m, depth=depth, zoom=zoom)
    fig.savefig('temp.png', bbox_inches='tight')
        # import pdb; pdb.set_trace()
