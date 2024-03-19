import math
import os
import logging
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as num
import gc
from pyrocko.guts import load
from pyrocko.plot.automap import Map
from pyrocko import util
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_allgains(self, results_all, stats_list, directory, fn):
    n_stats = results_all.shape[1]
    n_plotrows = math.ceil(n_stats/50.)+1
    ax = ['ax'+str(nrow) for nrow in range(n_plotrows)]
    ax = tuple(ax)
    fig, ax = plt.subplots(nrows=n_plotrows, figsize=(20, n_plotrows*5))

    maxy = num.nanmax(results_all)
    miny = num.nanmin(results_all)

    times = [util.time_to_str(ev.time)[:10] for ev in self.events]

    cmap = plt.get_cmap('viridis')
    indices = num.linspace(0, cmap.N, len(self.events))
    my_colors = [cmap(int(i)) for i in indices]

    bounds = num.linspace(0, len(times), len(times)+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    # sort result_all (axis = 2) und stat_list nach stat_list
    nos_list = [i for i in range(len(stats_list))]
    ind_stats_list_sorted = [x for _, x in sorted(zip(list(stats_list), nos_list))]
    stats_list = sorted(list(stats_list))
    results_all_sorted = num.empty((results_all.shape))

    for i, j in zip(ind_stats_list_sorted, range(len(ind_stats_list_sorted))):
        results_all_sorted[:, j] = results_all[:, i]

    for i, (row, axes) in enumerate(zip(range(n_plotrows), list(ax))):
        # if i != n_plotrows:
        start = row*50
        if start+50 < n_stats:
            stop = start+50
        else:
            stop = n_stats
        for i_st in (range(start, stop)):
            yval = results_all_sorted[:, i_st]
            x = i_st - row*50
            xval = [x for val in range(len(yval))]

            # make the scatter
            axes.scatter(xval, yval, c=my_colors, s=60, marker='o', cmap=cmap,
                         norm=norm, edgecolor='none')

            # axes.set_xlabel('Station', fontsize=14)
            axes.set_yscale('log')
            if len(self.method) == 2:
                axes.set_ylabel('Gain relative to %s %s' % (str(self.method[1][0]),
                                str(self.method[1][1])), fontsize=14)
            elif self.method == 'scale_one':
                axes.set_ylabel('Absolute gain')
            else:
                axes.set_ylabel('Relative gain')
            axes.set_ylim(miny, maxy)

            stats = [str(st[0])+' '+str(st[1])
                     for st in list(stats_list)[start:stop]]
            if self.method == 'syn':
                stats = stats_list[start:stop]
            axes.set_xticks(num.arange(len(stats)))
            axes.set_xticklabels(stats, rotation=60)

        # if i == n_plotrows:
    # ax2 = fig.add_axes([0.3, 0.1, 0.4, 0.01])
    # plt.axis('off')
    ax[-1].set_xticks([])
    ax[-1].set_yticks([])
    ax[-1].spines['top'].set_visible(False)
    ax[-1].spines['right'].set_visible(False)
    ax[-1].spines['bottom'].set_visible(False)
    ax[-1].spines['left'].set_visible(False)

    divider = make_axes_locatable(ax[-1])
    cax = divider.append_axes('top', size='5%', pad=0.1)
    # plt.subplots_adjust(bottom=0.9)
    # cax = plt.axes([0.3, 0.1, 0.4, 0.01])
    cbar = mpl.colorbar.ColorbarBase(ax=cax, cmap=cmap, norm=norm,
                                     spacing='uniform',
                                     ticks=bounds,
                                     boundaries=bounds,
                                     orientation='horizontal')
    # plt.colorbar(cax=cax)
    try:
        cbar.ax.set_xticklabels(times, rotation=60)
    except:
        # work around for newer matplotlib versions...
        times.append('')
        cbar.ax.set_xticklabels(times, rotation=60)
    # ax[-1].axis('off')

    plt.tight_layout()
    plt.savefig(os.path.join(directory, fn)[:-4]+'_all_gains.png')


def plot_median_gain_map_from_file(ns,
                                   st_lats,
                                   st_lons,
                                   pl_options,
                                   pl_topo,
                                   gains_file,
                                   directory,
                                   comp,
                                   mapsize,
                                   outformat):
    """
    Plot map with mean relative gain factors

    """
    gmtconf = dict(
                   MAP_TICK_PEN_PRIMARY='1.25p',
                   MAP_TICK_PEN_SECONDARY='1.25p',
                   MAP_TICK_LENGTH_PRIMARY='0.2c',
                   MAP_TICK_LENGTH_SECONDARY='0.6c',
                   FONT_ANNOT_PRIMARY='20p,1,black',
                   FONT_LABEL='20p,1,black',
                   PS_CHAR_ENCODING='ISOLatin1+',
                   MAP_FRAME_TYPE='fancy',
                   FORMAT_GEO_MAP='D',
                   PS_PAGE_ORIENTATION='portrait',
                   MAP_GRID_PEN_PRIMARY='thinnest,0/50/0',
                   MAP_ANNOT_OBLIQUE='6')

    gains_fromfile = load(filename=os.path.join(directory, gains_file))
    try:
        ns_rel = gains_fromfile.ref_stats
        # print(ns_rel)
    except Exception:
        ns_rel = None

    gains_no_nan = []
    lat_no_nan = []
    lon_no_nan = []
    # stats_no_nan = []

    # print(gains_fromfile.trace_gains_median)

    for i_ns, ns_now in enumerate(ns):
        for l in ['00', '', '01', '10', '60']:
            try:
                nslc = '%s.%s.%s.%s' % (ns_now[0], ns_now[1], l, comp)
                g = gains_fromfile.trace_gains_median[nslc]

                if g > 0.0 or g < 0.0:  # g < 5 and g > 0.2:  # g < 10.0 and g > 0.1: #g > 0.0 or g < 0.0:
                    gains_no_nan.append(g)
                    # stats_no_nan.append(ns_now[1])
                    lat_no_nan.append(st_lats[i_ns])
                    lon_no_nan.append(st_lons[i_ns])
            except KeyError:
                continue

    miny = min(gains_no_nan)
    maxy = max(gains_no_nan)
    # print(gains_no_nan)
    gains_fromfile = None
    gc.collect()
    gains_no_nan = list(num.log10(gains_no_nan))
    # print(gains_no_nan)
    m = Map(
        lat=pl_options[0],
        lon=pl_options[1],
        radius=pl_options[2],
        width=mapsize[0],
        height=mapsize[1],
        show_grid=False,
        show_topo=pl_topo,
        # topo_cpt_dry='/home/gesap/Documents/CETperceptual_GMT/CET-L2.cpt',#'/usr/local/share/cpt/gray.cpt',
        color_dry=(143, 188, 143),  # grey
        illuminate=True,
        illuminate_factor_ocean=0.15,
        # illuminate_factor_land = 0.2,
        show_rivers=True,
        show_plates=False,
        gmt_config=gmtconf)

    # Draw some larger cities covered by the map area
    # m.draw_cities()

    # Draw max. amplitudes at station locations as colored circles
    cptfile = 'tempfile.cpt'
    miny = min(gains_no_nan)
    maxy = max(gains_no_nan)

    # plus/minus 10 % for nicer scale
    miny = miny-0.1*abs(miny)
    maxy = maxy+0.1*abs(maxy)

    try:
        m.gmt.makecpt(
            C=pl_options[3],#'split',#'polar',  # '/home/gesap/Documents/CETperceptual_GMT/CET-D4.cpt',#'split',#'polar',
            T='%f/%f' % (miny, maxy),  # (miny, maxy), # (-1,1),#(-0.7, 0.7), (-20, 20)
            Q=True,
            out_filename=cptfile, suppress_defaults=True)
    except:
        try:
            m.gmt.makecpt(
                C='split',
                T='%f/%f' % (miny, maxy),
                Q=True,
                out_filename=cptfile, suppress_defaults=True)
            logging.warning('Could not find gmt cptfile, using split instead.')
        except:
            logging.error('Could not make gmt cpt file for map.')

    # if ns_rel:
    #     text = 'Gains relative to station %s, colorbar values 10^x.'\
    #            % (ns_rel)
    # else:
    #     text = 'colorbar values 10^x'
    # m.gmt.psbasemap(B='+t"%s"' % text, *m.jxyr)

    m.gmt.psxy(in_columns=(lon_no_nan, lat_no_nan, gains_no_nan),
               C=cptfile, G='blue', S='t24p',
               *m.jxyr)

    # add a colorbar
    B_opt_psscale = 'xaf'
    m.gmt.psscale(
                B=B_opt_psscale+'+l log10(Md(A@-i,j@-/A@-ref,j@-))',
                D='x9c/6c+w12c/0.5c+jTC+h',
                C=cptfile)

    # add station labels

    # for i in range(len(stats_no_nan)):
    #    m.add_label(lat_no_nan[i], lon_no_nan[i], stats_no_nan[i])
    fn = os.path.join(directory, '%s%s_map_log.%s' % (gains_file[0:12], comp, outformat))
    m.save(fn)
    logging.info('saved file %s' % fn)
