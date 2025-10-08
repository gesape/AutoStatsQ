import glob, os, logging
from pyrocko import guts 


logger = logging.getLogger('REPORT-ORIENT')


def include_map(orient_result_dir, ccmin):

    orient_map = os.path.join(orient_result_dir, 'map_orient_%s.png' % ccmin)
    print(orient_map)

    if not glob.glob(orient_map):

        logger.info('Error: orient map not found.')
        error = ''' \t \t \t \t<section data-transition="fade">\n
                    \t \t \t \t \t <h4>Overview map of orientation results:</h4>\n
                    \t \t \t \t \t<p style="font-size:%s">Map not found.</p>\n
                    \t \t \t \t</section>
                    ''' % ('50%'.format())

        return error

    else:

        include_this = ''' \t \t \t \t<section data-transition="fade">\n
        \t \t \t \t \t <h4>Overview map of orientation results (ccmin=%s):</h4>\n
        \t \t \t \t \t<img style="margin:0; margin-bottom:-0.6em; margin-right:0.3em; height: 12em;"
                                src="../results/orient/map_orient_%s.png" />\n
        \t \t \t \t \t<p style="font-size:%s">Orientation of the N component obtained from median correction angles. Small black arrows indicate station for which less than 5 events could be used.</p>\n
        \t \t \t \t</section>
        ''' % (ccmin,ccmin,'50%'.format())

        return include_this

          

def include_table_correction_angles(orient_result_dir, ccmin):
    try:
        rota_stats = guts.load(filename='%s/CorrectionAngles_cc%s.yaml' % (orient_result_dir, ccmin))
    except FileNotFoundError:
        logger.info('Error: Correction angle table not found.')
        error = ''' \t \t \t \t<section data-transition="fade">\n
                    \t \t \t \t \t <h4>Orientation test results - statistics</h4>\n
                    \t \t \t \t \t<p style="font-size:%s">File not found.</p>\n
                    \t \t \t \t</section>
                    ''' % ('50%'.format())

        return error

    medians = rota_stats.CorrectAngl_perStat_median
    means = rota_stats.CorrectAngl_perStat_mean
    stds = rota_stats.CorrectAngl_perStat_stdd
    nevs = rota_stats.n_events

    stats_list = []

    for key in medians.keys():

        med = round(medians[key],1)
        mean = round(means[key],1)
        std = round(stds[key],1)
        nev = nevs[key]

        if nev < 5:
            single_stat_in = '''
                    <tr>\n
                        <td style="color:#c8c8c8">%s</td>\n
                        <td style="color:#c8c8c8">%s</td>\n
                        <td style="color:#c8c8c8">%s</td>\n
                        <td style="color:#c8c8c8">%s</td>\n
                        <td style="color:#c8c8c8">%s</td>\n
                    </tr>\n''' % (key, med, mean, std, nev)

        elif abs(med) > 10:
            single_stat_in = '''
                    <tr>\n
                        <td style="color:#e61717">%s</td>\n
                        <td style="color:#e61717">%s</td>\n
                        <td style="color:#e61717">%s</td>\n
                        <td style="color:#e61717">%s</td>\n
                        <td style="color:#e61717">%s</td>\n
                    </tr>\n''' % (key, med, mean, std, nev)

        else:

            single_stat_in = '''
                    <tr>\n
                        <td>%s</td>\n
                        <td>%s</td>\n
                        <td>%s</td>\n
                        <td>%s</td>\n
                        <td>%s</td>\n
                    </tr>\n''' % (key, med, mean, std, nev)

        
        stats_list.append(single_stat_in)

    stats = ''
    for s in stats_list:
        stats += s
    #print(stats)


    #sys.exit()
    string1 = '''
    <section>
        <h4>Orientation test results - statistics</h4>
        <p style="font-size:50%">Grey font color for stations with less than 5 stable results, red font color for misorientations larger 10°.</p>\n
        <div style="overflow:scroll; height:500px; font-size: 40%;">
        <table>
            <thead><tr>
                <th>Station</th>
                <th>Median</th>
                <th>Mean</th>
                <th>Std dev</th>
                <th>n events</th>
            </tr></thead>
            <tbody>''' #% ('40%'.format())
    #print(string1)
    string2 = stats
    string3 = '''
            </tbody>
        </table>
        </div>
    </section>
    '''
    string = string1+string2+string3
    return string



def add_all_stats_figures_timing(orient_result_dir, ccmin):

    timing_plots = glob.glob(os.path.join(orient_result_dir, '*_overtime_%s.png' % ccmin))

    if len(timing_plots) <1:

        logger.info('Error: orient over time plots not found.')
        error = ''' \t \t \t \t<section data-transition="fade">\n
                    \t \t \t \t \t <h4>Median orientation angle over time</h4>\n
                    \t \t \t \t \t<p style="font-size:%s">Figure not found.</p>\n
                    \t \t \t \t</section>
                    ''' % ('50%'.format())

        return error

    else:

        str1 = '''
                <section>
                    <h4>Median orientation angle over time</h4>
                    <div style="overflow:scroll; height:500px; font-size: 30%;">\n'''
        str2 = ''
        for ot in timing_plots:
            pl = ot.split('/')[-1]
            str2 += '''<img style="margin:0; margin-bottom:-0.6em; margin-right:0.3em; height: 12em;"
                                src="../results/orient/%s" />\n''' % (pl)
        str3 = '''
                </section>'''

    return str1+str2+str3



def add_all_stats_figures_baz(orient_result_dir, ccmin):

    baz_plots = glob.glob(os.path.join(orient_result_dir, '*_baz_%s.png' % ccmin))

    if len(baz_plots) <1:

        logger.info('Error: orient over backazimuth plots not found.')
        error = ''' \t \t \t \t<section data-transition="fade">\n
                    \t \t \t \t \t <h4>Median orientation angle over time</h4>\n
                    \t \t \t \t \t<p style="font-size:%s">Figure not found.</p>\n
                    \t \t \t \t</section>
                    ''' % ('50%'.format())

        return error

    else:

        str1 = '''
                <section>
                    <h4>Median orientation angle vs. event backazimuth</h4>
                    <div style="overflow:scroll; height:500px; font-size: 30%;">\n'''
        str2 = ''
        for ot in baz_plots:
            pl = ot.split('/')[-1]
            str2 += '''<img style="margin:0; margin-bottom:-0.6em; margin-right:0.3em; height: 12em;"
                                src="../results/orient/%s" />\n''' % (pl)
        str3 = '''
                </section>'''

    return str1+str2+str3



def add_all_stats_figures_ccs(orient_result_dir, ccmin):

    cc_plots = glob.glob(os.path.join(orient_result_dir, '*_distr_%s.*' % ccmin))

    if len(cc_plots) <1:
        logger.info('Error: orient cc plots not found.')
        error = ''' \t \t \t \t<section data-transition="fade">\n
                    \t \t \t \t \t <h4>Cross-correlation values over rotation angle, for each station and event</h4>\n
                    \t \t \t \t \t<p style="font-size:%s">Figure not found.</p>\n
                    \t \t \t \t</section>
                    ''' % ('50%'.format())

        return error

    else:
        str1 = '''
                <section>
                    <h4>Cross-correlation values over rotation angle, for each station and event</h4>
                    <div style="overflow:scroll; height:500px; font-size: 40%;">
                    <p>For each station, subplots for all events are shown. Empty plots indicate missing data. The cross-correlation value is plotted over correction angle. A correctly oriented sensor corresponds to a curve with a clear maximum at a correction angle of 0° and minima at -180/180°.</p>\n
    \n'''
        str2 = ''
        for ot in cc_plots:
            path = ot.split('/')
            plo = path[-1]
            pnet,stat,loc,pl,ccmin = path[-1].split('_')
            st = '%s.%s' % (pnet, stat)
            #str2 += '''<img style="margin:0; margin-bottom:-0.6em; margin-right:0.3em; height: 20em;"
            #                    src="%s" />\n''' % (ot)
            str2 +=         '''<p>%s</p>\n''' % st
            str2 +=         '''<embed src="../results/orient/%s" />\n''' % plo
        str3 = '''
                </section>'''

        return str1+str2+str3


def add_single_station_results(orient_result_dir, ccmin):

    try:
        rota = guts.load(filename='%s/AllCorrectionAngles_cc%s.yaml' % (orient_result_dir, ccmin)).dict_stats_all
    except FileNotFoundError:
        logger.info('Error: Correction angle single station table not found.')
        error = ''' \t \t \t \t<section data-transition="fade">\n
                    \t \t \t \t \t <h4>Orientation test - single station results</h4>\n
                    \t \t \t \t \t<p style="font-size:%s">File not found.</p>\n
                    \t \t \t \t</section>
                    ''' % ('50%'.format())

        return error

    
    string1 = '''
            <section>
                <h4>Orientation test - single station results</h4>
                <p style="font-size:50%">Tip: Use str f to search for results of a single station.</p>\n
                <div style="overflow:scroll; height:500px; font-size: 30%;">\n'''

    string2 = ''
    for dtdict in rota:
        if dtdict.ev_rota:
            stat = '%s.%s' % (dtdict.station[0],dtdict.station[1])
            string2 += '''
                    <p style="font-size:120%">''' + '''%s</p>\n''' % stat
            string2 += '''
                    <table>\n
                    <thead><tr>\n
                    <th>Event date</th>\n
                    <th>Correction angle</th>\n
                    </tr></thead>\n
                    <tbody>\n'''
            
            string3 = ''

            items = dtdict.ev_rota.items()
            for evt,a in sorted(items):

                string3+='''
                    <tr>\n''' +\
                    '''<td>%s</td>\n''' % evt 
                string3 += '''
                    <td>%s</td>\n''' % str(a)

            string2 += string3 +'''
                </tr>\n
                </tbody>\n
                </table>\n
                
                ''' 



    string4 = '''
    </div>\n        
    </section>
    '''

    string = string1+string2+string4

    return string

