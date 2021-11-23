import sys, glob, os, logging
from pyrocko import guts, model, util

logger = logging.getLogger('REPORT-CATALOG')


def include_maps_subsets(tmin, tmax, catalog_result_dir):
    deep_fig = glob.glob(os.path.join(catalog_result_dir, 'catalog_*deep_subset.png'))
    shallow_fig = glob.glob(os.path.join(catalog_result_dir, 'catalog_*shallow_subset.png'))
    try:
        deep_fig = deep_fig[0]
        shallow_fig = shallow_fig[0]
    except:
        logger.info('ERROR - catalog subset maps not found.')
        error = '''
        <section>
            <h4>Catalog subset maps</h4>
            <p style="font-size:100%">Files not found.</p>\n
        </section>
        '''
        return error

    deep_fig = os.path.join('../results/catalog', deep_fig.split('/')[-1])
    shallow_fig = os.path.join('../results/catalog', shallow_fig.split('/')[-1])

    include_this = ''' \t \t \t \t<section data-transition="fade">\n
    \t \t \t \t \t <h4>Subsets of teleseismic events from GCMT catalog:</h4>\n
    \t \t \t \t \t <p style="font-size:%s">%s to %s:</p>\n
    \t \t \t \t \t<img style="margin:0; margin-bottom:-0.6em; margin-right:0.3em; height: 8em;"
                            src="%s" />\n
    \t \t \t \t \t<img style="margin:0; margin-bottom:-0.6em; margin-right:0.3em; height: 8em;"
                            src="%s" />\n                            
    \t \t \t \t \t<p style="font-size:%s">Deep (left) and shallow (right) subsets of teleseismic events used for the AutoStatsQ tests.</p>\n
    \t \t \t \t</section>
    ''' % ('80%'.format(), tmin[:10], tmax[:10], deep_fig, shallow_fig, '40%'.format())

    return include_this


def include_event_tabels(tmin, tmax, Mmin, catalog_result_dir):
    deep_catf = glob.glob('%s/catalog_Mgr%s_deep.txt' % (catalog_result_dir, Mmin))
    shallow_catf = glob.glob('%s/catalog_Mgr%s_shallow.txt' % (catalog_result_dir, Mmin))

    try:
        deep_cat = model.load_events(deep_catf[0])
        shallow_cat = model.load_events(shallow_catf[0])
    except:
        error = '''
        <section>
            <h4>Catalog subset tables</h4>
            <p style="font-size:100%">Files not found.</p>\n
        </section>
        '''
        return error

    ev_list = []
    for cat in ['deep', 'shallow']:
        if cat == 'deep':
            events = deep_cat
        else:
            events = shallow_cat

        for ev in events:

            single_ev_in = '''
                    <tr>\n
                        <td>%s</td>\n
                        <td>%s</td>\n
                        <td>%s</td>\n
                        <td>%s</td>\n
                        <td>%s</td>\n
                        <td>%s</td>\n
                    </tr>\n''' % (cat, util.tts(ev.time), round(ev.lat,1), round(ev.lon,1), round(ev.magnitude,1), ev.region)

            
            ev_list.append(single_ev_in)

    evs = ''
    for s in ev_list:
        evs += s
    #print(stats)

    #sys.exit()
    string1 = '''
    <section>
        <p style="font-size:40%">Event catalogs</p>\n
        <div style="overflow:scroll; height:500px; font-size: 40%;">
        <table>
            <thead><tr>
                <th>Subset</th>
                <th>Date/Time</th>
                <th>Lat.</th>
                <th>Lon.</th>
                <th>M</th>
                <th>Region</th>
            </tr></thead>
            <tbody>''' #% ('40%'.format())
    #print(string1)
    string2 = evs
    string3 = '''
            </tbody>
        </table>
        </div>
    </section>
    '''
    string = string1+string2+string3
    return string


def include_statistic_plots(catalog_result_dir):
    deep_fig = glob.glob('%s/cat_hist_magn_dist_deep_subset.png' % catalog_result_dir)
    shallow_fig = glob.glob('%s/cat_hist_magn_dist_shallow_subset.png' % catalog_result_dir)
    try:
        deep_fig = deep_fig[0]
        shallow_fig = shallow_fig[0]
    except:
        logger.info('ERROR - catalog statistic plots not found.')
        error = '''
        <section>
            <h4>Catalog statistics plots</h4>
            <p style="font-size:100%">Files not found.</p>\n
        </section>
        '''
        return error

    deep_fig = os.path.join('../results/catalog', deep_fig.split('/')[-1])
    shallow_fig = os.path.join('../results/catalog', shallow_fig.split('/')[-1])

    include_this = ''' \t \t \t \t<section data-transition="fade">\n
    \t \t \t \t \t <p style="font-size:%s">Statistic plots of catalog subsets:</p>\n
    \t \t \t \t \t<img style="margin:0; margin-bottom:-0.6em; margin-right:0.3em; height: 8em;"
                            src="%s" />\n
    \t \t \t \t \t<img style="margin:0; margin-bottom:-0.6em; margin-right:0.3em; height: 8em;"
                            src="%s" />\n                            
    \t \t \t \t \t<p style="font-size:%s">Deep (left) and shallow (right) subsets of teleseismic events used for the AutoStatsQ tests.</p>\n
    \t \t \t \t</section>
    ''' % ('60%'.format(), deep_fig, shallow_fig, '40%'.format())

    return include_this

