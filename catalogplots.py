import numpy as num
from pyrocko import gmtpy
import matplotlib.pyplot as plt


def gmtplot_catalog_azimuthal(catalog, mid_point, dist, outfile):
    '''
    Plot events of catalog on map

    :param catalog: event catalog in pyrocko format
    :param mid_point: centre of map
    :param dist: max. distance in degrees
    '''
    gmt = gmtpy.GMT(config={'MAP_GRID_PEN_PRIMARY': '0.1p',
                            'MAP_FRAME_TYPE': 'fancy',
                            'FONT_ANNOT_PRIMARY': '8p,Helvetica,black',
                            'FONT_ANNOT_SECONDARY': '8p,Helvetica,black',
                            'FONT_LABEL': '8p,Helvetica,black'})
    gmt.psbasemap(
                  R='0/360/-90/0',
                  J='S0/-90/90/5i',
                  B='xa15')
    gmt.pscoast(
                R='g',
                J='E%s/%s/%s/5i' % (str(mid_point[1]),
                                    str(mid_point[0]),
                                    str(dist)),
                D='c',
                G='darkgrey')
    gmt.psxy(R='0/360/-90/0',
             J='E%s/%s/%s/5i' % (str(mid_point[1]),
                                 str(mid_point[0]),
                                 str(dist)),
             in_columns=([ev.lon for ev in catalog],
                         [ev.lat for ev in catalog]),
             G='red',
             S='a5p')
    gmt.psbasemap(
                  R='g',
                  J='E0/-90/%s/5i' % (str(dist)),
                  B='xg15')
    gmt.psbasemap(
                  R='g',
                  J='E0/-90/%s/5i' % (str(dist)),
                  B='yg15')
    gmt.save(outfile)
    print('saved ', outfile)


def plot_catalog_hist(catalog, dist_array, mean_wedges_mp, bins_hist, data_dir,
                      min_mag,
                      plot_wedges_vs_dist=True,
                      plot_wedges_vs_magn=True):
    '''
    Plot histogram showing number of earthquakes in each 15deg wedge
    '''
    mean_dists = num.mean(dist_array, axis=1)
    fig, ax1 = plt.subplots()
    ax1.hist(mean_wedges_mp, normed=False, bins=bins_hist, color='lightgrey')
    ax1.set_xlabel('15' + r'$^\circ$' + ' backazimuth bins')
    ax1.set_ylabel('Number of events in catalog')
    ax1.set_xlim((0, 24))
    ax1.set_ylim((0.1, 100))
    ax1.set_yscale('log', nonposy='clip')

    if plot_wedges_vs_dist is True:
        '''
        Plot distance along with frequency
        '''
        ax2 = ax1.twinx()
        ax2.plot(mean_wedges_mp+0.5, mean_dists/1000, 'ro')
        ax2.tick_params('y', colors='r')
        ax2.set_ylabel('Distance (km)', color='r')
        ax2.set_ylim((0, 30000))

    if plot_wedges_vs_magn is True:
        '''
        Plot magnitude along with frequency
        '''
        ax3 = ax1.twinx()
        if plot_wedges_vs_dist is True:
            ax3.spines['right'].set_position(('axes', 1.2))
        ax3.plot(mean_wedges_mp+0.5, [ev.magnitude for ev in catalog], 'bo')
        ax3.tick_params('y', colors='b')
        ax3.set_ylim(min_mag-0.2, 10)
        ax3.set_ylabel('Magnitude '+r'$(M_W)$', color='b')
    fig.tight_layout()
    plt.savefig(data_dir+'results/catalog/cat_hist_magn_dist.png')


def plot_distmagn(dist_array, catalog, data_dir):
    mean_dists = num.mean(dist_array, axis=1)
    fig, ax = plt.subplots()
    ax.plot(mean_dists/1000, [ev.magnitude for ev in catalog], 'ko')
    ax.set_xlabel('Distance '+r'$(km)$')
    ax.set_ylabel('Magnitude '+r'$(M_W)$')
    fig.tight_layout()
    plt.savefig(data_dir+'results/catalog/cat_dist_vs_magn.png')
