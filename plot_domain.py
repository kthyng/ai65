'''
Plots from Admiralty Inlet 65 meter simulation.
'''

import scipy.io
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from add_colormaps import add_colormaps
from skimage import color
from matplotlib import cm
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import MaxNLocator
import pdb
import cmocean.cm as cmo
import cmocean

mpl.rcParams.update({'font.size': 16})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'

def comps():
    '''read in comparison points'''

    import scipy.io
    base = 'savedoutput/'
    data = dict()
    d = scipy.io.loadmat(base + 'zeta_adcpzetasites.mat')
    data['zeta'] = dict(); data['vel'] = dict(); data['ctd'] = dict()
    data['zeta']['lon'] = d['coords']['xm'][0,0][0,:]; data['zeta']['lat'] = d['coords']['ym'][0,0][0,:]
    d = scipy.io.loadmat(base + 'uADCPvel.mat')
    data['vel']['lon'] = d['datacoords']['xm'][0,0][0,0,:]; data['vel']['lat'] = d['datacoords']['ym'][0,0][0,0,:]
    d = scipy.io.loadmat(base + 'PTH005_2006_9.mat')
    lon0 = d['PTH']['lon'][0][0][0][0]; lat0 = d['PTH']['lat'][0][0][0][0]
    d = scipy.io.loadmat(base + 'ADM001_2006_9.mat')
    lon1 = d['ADM']['lon'][0][0][0][0]; lat1 = d['ADM']['lat'][0][0][0][0]
    data['ctd']['lon'] = np.array([[lon0, lon1]]); data['ctd']['lat'] = np.array([[lat0, lat1]])

    # data markers
    data['zeta']['plot'] = {'color': '#FF5733', 'marker': 'o', 'markersize': 8, 'linewidth': 0}
    data['vel']['plot'] = {'color': 'k', 'marker': 's', 'markersize': 10,
                           'markerfacecolor': 'None', 'markeredgewidth': 0.8, 'linewidth': 0}
    data['ctd']['plot'] = {'color': '#8E44AD', 'marker': '*', 'markersize': 15, 'linewidth': 0}

    # d = scipy.io.loadmat(base + 'ha_zetaadcp.mat')
    # cons = d['ha']['name'][0][0]  # tidal constituents
    # freqs = d['ha']['freq'][0][0]  # tidal frequencies
    # # tidal constituents: constituents by 4:
    # tidecons = [d['ha']['tidecon'][0][0][i] for i in range(cons.size)]
    #
    # zeta = dict()
    # for con, freq, tidecon in zip(cons, freqs, tidecons):
    #     con = con.strip()  # remove trailing spaces
    #     zeta[con] = dict()
    #     zeta[con]['f'] = freq
    #     zeta[con]['tidecons'] = tidecon

    return data


def add_alpha(cmap, alpha):
    '''Add alpha to colormap. Return adjusted colormap.'''

    # Read in rgb of colormap
    rgb = cmocean.tools.print_colormaps([cmap], savefiles=False)[0]
    # make array of alpha values
    alphas = alpha*np.ones((rgb.shape[0],1))
    # set up new rgb
    rgbnew = np.hstack((rgb, alphas))
    # create new colormap object
    cmap = cmocean.tools.cmap(rgbnew)
    return cmap


def ai(grid, zoom='out'):
    '''
    Plot bathymetry of just Admiralty Inlet

    Inputs:
     grid       tracpy grid dictionary
     zoom       'out'= zoomed out view, 'in'= zoomed in view
    '''

    # x and y limits for this plot
    if zoom=='out': # zoomed out
        lonlims = [-122.8, -122.55]
        latlims = [47.9665, 48.227]
        # levels to plot
        levs = np.arange(-200, 20, 20)
    elif zoom=='in': # zoomed in
        lonlims = [-122.71, -122.64]
        latlims = [48.135, 48.18]
        # levels to plot
        levs = np.arange(-60, 5, 5)

    cmap = 'Blues_r'
    xlims, ylims = grid.proj(lonlims, latlims)

    # Make plot
    fig = plt.figure(figsize=(14,12))
    ax = fig.add_subplot(111)
    grid.proj.drawcoastlines(ax=ax)
    grid.proj.fillcontinents('darkgreen', ax=ax)
    mappable = ax.contourf(grid.x_rho, grid.y_rho, -grid.h, cmap=cmap, levels=levs, extend='min')
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    cb = fig.colorbar(mappable, shrink=0.75, pad=0.025)
    cb.set_label('Height/depth [m]')
    fig.show()

    # # Save figure
    # fig.savefig('figures/ai_transect.png')


# def ps(dosave=True, fname='figures/domains.png', lont=None, latt=None, ht=None, dd=None):
#     '''
#     Plot Bathymetry of Puget Sound, Admiralty Inlet, and Admiralty Head
#
#     Inputs:
#      dosave     Save figure
#      fname      File name for figure
#      lont, latt transect points to plot if they are input
#      ht         Depth along transect, if input
#      dd         Distance in meters along transect
#     '''

if __name__ == "__main__":  # Puget sound plot with inlays
    # download bathymetry, which can be found at: http://figshare.com/preview/_preview/1165560 (27.3MB)
    fname = 'figures/domains'
    dosave = True
    doAH = False
    lont=None; latt=None; ht=None; dd=None

    # Read in bathymetry
    mat = scipy.io.loadmat('../form-drag/cascadia_gridded.mat')

    # x and y limits for these plots
    lonlimsPS = np.array([-124.21, -122.15]) #-123.21, -122.15])
    latlimsPS = np.array([47.02, 48.82])
    lonlimsAI = np.array([-122.8, -122.535])
    latlimsAI = np.array([47.9665, 48.225])
    lonlimsAH = np.array([-122.72, -122.64])
    latlimsAH = np.array([48.12, 48.18])

    # Functionality copied from https://github.com/clawpack/geoclaw/blob/master/src/python/geoclaw/topotools.py#L873
    land_cmap = cmo.speed_r
    land_cmap = add_alpha(land_cmap, 0.7)
    sea_cmap = cmo.deep_r
    cmapPS, norm = add_colormaps((land_cmap, sea_cmap), data_limits=[-375,2500], data_break=0.0)
    cmapAI = cmo.deep_r
    cmapAH = cmo.deep_r

    cland = '#4C6351'

    # levels to plot
    levsPS = np.concatenate((np.arange(-375, 0, 25), np.arange(0,3000,500)))
    levsAI = np.arange(-200, 20, 20)
    levsAH = np.arange(-120, 15, 15)

    # use basemap
    basemapPS = Basemap(llcrnrlon=lonlimsPS[0], llcrnrlat=latlimsPS[0],
                    urcrnrlon=lonlimsPS[1], urcrnrlat=latlimsPS[1],
                    lat_0=latlimsPS.mean(), lon_0=lonlimsPS.mean(),
                    projection='lcc', resolution='f',
                    area_thresh=0.)
    xPS, yPS = basemapPS(mat['lon_topo'], mat['lat_topo'])
    xlimsAI, ylimsAI = basemapPS(lonlimsAI, latlimsAI)
    xlimsAH, ylimsAH = basemapPS(lonlimsAH, latlimsAH)

    # data points to plot on Admiralty Inlet map, get projected coords
    data = comps()
    for key in data.keys():
        data[key]['x'], data[key]['y'] = basemapPS(data[key]['lon'], data[key]['lat'])

    # Make Puget Sound plot
    fig = plt.figure(figsize=(16,16))
    axPS = fig.add_subplot(111)
    # basemapPS.drawcoastlines(ax=axPS)
    mappablePS = axPS.contourf(xPS, yPS, mat['z_topo'], cmap=cmapPS, levels=levsPS, zorder=2, norm=norm)
    locator = MaxNLocator(6) # if you want no more than 10 contours
    locator.create_dummy_axis()
    locator.set_bounds(lonlimsPS[0], lonlimsPS[1])
    pars = locator()
    locator = MaxNLocator(6) # if you want no more than 10 contours
    locator.create_dummy_axis()
    locator.set_bounds(latlimsPS[0], latlimsPS[1])
    mers = locator()
    basemapPS.drawparallels(mers, dashes=(1, 1), linewidth=0.15, labels=[1,0,0,0], ax=axPS)#, zorder=3)
    basemapPS.drawmeridians(pars, dashes=(1, 1), linewidth=0.15, labels=[0,0,0,1], ax=axPS)#, zorder=3)
    cbPS = fig.colorbar(mappablePS, pad=0.015, aspect=35)
    cbPS.set_label('Height/depth [m]')
    # Label
    axPS.text(0.8, 0.025, 'Puget Sound', transform=axPS.transAxes, color='0.15')

    # Inset magnified plot of Admiralty Inlet
    locAI = 3
    axAI = zoomed_inset_axes(axPS, 4, loc=locAI)
    basemapPS.drawcoastlines(ax=axAI)
    basemapPS.fillcontinents(cland, ax=axAI)
    mappableAI = axAI.contourf(xPS, yPS, mat['z_topo'], cmap=cmapAI, levels=levsAI)
    # add data points
    for key, value in zip(data.keys(), data.values()):
        xtemp, ytemp = data[key]['x'], data[key]['y']
        axAI.plot(xtemp, ytemp, **data[key]['plot'])
    axAI.set_xlim(xlimsAI)
    axAI.set_ylim(ylimsAI)
    # Inlaid colorbar
    if locAI == 1:
        caxAI = fig.add_axes([0.581, 0.665, 0.011, 0.1])
    elif locAI == 3:
        caxAI = fig.add_axes([0.18, 0.12, 0.011, 0.1])
    cbAI = plt.colorbar(mappableAI, cax=caxAI, orientation='vertical')
    cbAI.ax.tick_params(labelsize=12)
    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    mark_inset(axPS, axAI, loc1=2, loc2=4, fc="none", ec="0.3", lw=1.5, zorder=5)
    # Label
    axAI.text(0.41, 0.83, 'Admiralty\n      Inlet', transform=axAI.transAxes, color='0.15', fontsize=16)

    # Inset magnified plot of Admiralty Head
    if doAH:
        axAH = zoomed_inset_axes(axPS, 9, loc=1)
        basemapPS.drawcoastlines(ax=axAH)
        basemapPS.fillcontinents(cland, ax=axAH)
        mappableAH = axAH.contourf(xPS, yPS, mat['z_topo'], cmap=cmapAH, levels=levsAH)
        axAH.set_xlim(xlimsAH)
        axAH.set_ylim(ylimsAH)

        if plt.is_numlike(lont):
                # add points if you have some
                xt, yt = basemapPS(lont, latt)
                axAH.plot(xt, yt, 'k', lw=3)

            # Inlaid colorbar
        caxAH = fig.add_axes([0.398, 0.116, 0.012, 0.15])
        cbAH = plt.colorbar(mappableAH, cax=caxAH, orientation='vertical')
        cbAH.ax.tick_params(labelsize=12)
        # draw a bbox of the region of the inset axes in the parent axes and
        # connecting lines between the bbox and the inset axes area
        mark_inset(axPS, axAH, loc1=2, loc2=4, fc="none", ec="0.3", lw=1.5, zorder=5)
        # Label
        axAH.text(0.47, 0.92, 'Admiralty Head', transform=axAH.transAxes, color='0.15', fontsize=16)

    # pdb.set_trace()

    if plt.is_numlike(lont):
        # Add axes to plot transect depths
        axdepths = fig.add_axes([0.28, 0.39, 0.14, 0.075], zorder=11)
        axdepths.plot((np.arange(lont.size)*dd)/1000., -ht, '0.2', lw=2, zorder=12)
        axdepths.tick_params(axis='both', colors='0.1', top='off', right='off', width=2, length=4, labelsize=12, labelcolor='0.1')
        axdepths.spines['bottom'].set_color('none')
        axdepths.spines['top'].set_color('none')
        axdepths.spines['left'].set_color('none')
        axdepths.spines['right'].set_color('none')
        axdepths.set_xlabel('Distance along transect [km]', fontsize=14, color='0.1')
        axdepths.set_ylabel('Transect depth [m]', fontsize=14, color='0.1')
        axdepths.patch.set_alpha(0.0) # make bg transparent
        fig.show()


    # Save figure
    if dosave:
        fig.savefig(fname + '.png', bbox_inches='tight', dpi=120)
        fig.savefig(fname + '_highres.png', bbox_inches='tight', dpi=300)
    fig.show()
