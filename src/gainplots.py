import math
import statistics
import matplotlib.pyplot as plt
import numpy as num
import gc
from pyrocko.guts import load
from pyrocko.plot.automap import Map


def plot_allgains(self, results_all, stats_list, directory, fn):
    n_stats = results_all.shape[1]
    n_plotrows = math.ceil(n_stats/50.)
    ax = ['ax'+str(nrow) for nrow in range(n_plotrows)]
    ax = tuple(ax)
    fig, ax = plt.subplots(nrows=n_plotrows, figsize=(20, n_plotrows*5))

    maxy = num.nanmax(results_all)
    miny = num.nanmin(results_all)

    for row, axes in zip(range(n_plotrows), list(ax)):
        start = row*50
        if start+50 < n_stats:
            stop = start+50
        else:
            stop = n_stats
        for i_st in (range(start, stop)):
            yval = results_all[:, i_st]
            x = i_st - row*50
            xval = [x for val in range(len(yval))]
            axes.plot(xval, yval, '.')
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
    plt.tight_layout()
    plt.savefig(directory+fn[:-4]+'_all_gains.png')


def plot_mean_gain_map_from_file(ns,
                                 st_lats,
                                 st_lons,
                                 pl_options,
                                 gains_file,
                                 directory,
                                 comp):
    '''
    Plot map with mean relative gain factors

    '''
    gains_fromfile = load(filename=directory+gains_file)
    try:
        ns_rel = gains_fromfile.ref_stats
        print(ns_rel)
    except:
        ns_rel = None

    gains_no_nan = []
    lat_no_nan = []
    lon_no_nan = []
    # stats_no_nan = []

    for i_ns, ns_now in enumerate(ns):
        try:
            nslc = '%s.%s..%s' % (ns_now[0], ns_now[1], comp)
            gains_no_nan.append(gains_fromfile.trace_gains[nslc])
            # stats_no_nan.append(ns_now[1])
            lat_no_nan.append(st_lats[i_ns])
            lon_no_nan.append(st_lons[i_ns])
        except KeyError:
            continue
    gains_fromfile = None
    gc.collect()
    gains_no_nan = list(num.log10(gains_no_nan))

    m = Map(
        lat=pl_options[0],
        lon=pl_options[1],
        radius=pl_options[2],
        width=30., height=30.,
        show_grid=False,
        show_topo=False,
        # show_topo_scale=True,
        # color_dry=(238, 236, 230),
        # topo_cpt_wet='light_sea_uniform',
        # topo_cpt_dry='light_land_uniform',
        illuminate=True,
        illuminate_factor_ocean=0.15,
        show_rivers=True,
        show_plates=True)

    # Draw some larger cities covered by the map area
    m.draw_cities()

    # Draw max. amplitudes at station locations as colored circles
    cptfile = 'tempfile.cpt'
    miny = min(gains_no_nan)
    maxy = max(gains_no_nan)

    # plus/minus 10 % for nicer scale
    miny = miny-0.1*miny
    maxy = maxy+0.1*maxy
    m.gmt.makecpt(
                C='rainbow',
                T='%g/%g' % (miny, maxy),
                Q=True,
                out_filename=cptfile, suppress_defaults=True)

    if ns_rel:
        text = 'Gains relative to station %s, colorbar values 10^x.'\
               % (ns_rel)
    else:
        text = 'colorbar values 10^x'
    m.gmt.psbasemap(B='+t"%s"' % text, *m.jxyr)

    m.gmt.psxy(in_columns=(lon_no_nan, lat_no_nan, gains_no_nan),
               C=cptfile, G='blue', S='t10p',
               *m.jxyr)

    # add a colorbar
    B_opt_psscale = 'xaf'
    m.gmt.psscale(
                B=B_opt_psscale,
                D='x9c/4c+w12c/0.5c+jTC+h',
                C=cptfile)

    # add station labels

    # for i in range(len(stations)):
    #    m.add_label(lats[i], lons[i], stations[i])
    fn = directory+gains_file[:-4]+'_map_log.png'
    m.save(fn)
    print('saved file', fn)
