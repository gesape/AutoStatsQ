import glob, os, logging
from pyrocko import guts

logger = logging.getLogger('REPORT-GAINS')


def add_gain_map(comp, gain_result_dir):
    gain_map = os.path.join(gain_result_dir, 'gains_median%s_map_log.png' % (comp))

    if not glob.glob(gain_map):

        logger.info('Error: gain map not found for %s' % comp)
        error = ''' \t \t \t \t<section data-transition="fade">\n
                    \t \t \t \t \t <h4>Overview map of amplitude gain comparison results</h4>\n
                    \t \t \t \t \t <h4>%s component</h4>\n
                    \t \t \t \t \t<p style="font-size:%s">Map not found.</p>\n
                    \t \t \t \t</section>
                    ''' % (comp, '50%'.format())

        return error

    else:
        gain_map = os.path.join('../results/gains', gain_map.split('/')[-1])
        include_this = ''' \t \t \t \t<section data-transition="fade">\n
        \t \t \t \t \t <h4>Overview map of amplitude gain comparison results</h4>\n
        \t \t \t \t \t <h4>%s component</h4>\n
        \t \t \t \t \t<img style="margin:0; margin-bottom:0.3em; margin-right:0.3em; height: 9em;"
                                src="%s" />\n
        \t \t \t \t \t<p style="font-size:%s">Relative gains - median of ratios between max. P amplitude at signle stations compared to reference station or theoretical amplitude.</p>\n
        \t \t \t \t</section>
        ''' % (comp, gain_map, '50%'.format())

        return include_this


def add_gain_median_mean(comp, method, gain_result_dir):


    try:
        gain_stats = guts.load(filename='%s/gains_median_and_mean%s.txt' 
                            % (gain_result_dir, comp))
    except FileNotFoundError:
        logger.info('ERROR - Gain stats file not found.')
        error = '''
        <section>
            <h4>Amplitude gain ratio test results - statistics</h4>
            <p style="font-size:100%">File not found.</p>\n
        </section>
        '''
        return error

    else:
        medians = gain_stats.trace_gains_median
        means = gain_stats.trace_gains_mean
        stds = gain_stats.trace_gains_stdev
        nevs = gain_stats.n_ev_used

        if len(method) == 2:
            ref = '%s.%s' % (gain_stats.ref_stats.split(' ')[0], gain_stats.ref_stats.split(' ')[1])
        else:
            ref = ''
            method = method[0]

        stats_list = []

        for key in medians.keys():

            med = round(medians[key],1)
            mean = round(means[key],1)
            std = round(stds[key],1)
            nev = nevs[key]
            st = '%s.%s' % (key.split('.')[0], key.split('.')[1])

            if ref == st:  # here
                single_stat_in = '''
                        <tr>\n
                            <td style="color:#277545"><b>%s</b></td>\n
                            <td style="color:#277545"><b>%s</b></td>\n
                            <td style="color:#277545"><b>%s</b></td>\n
                            <td style="color:#277545"><b>%s</b></td>\n
                            <td style="color:#277545"><b>%s</b></td>\n
                        </tr>\n''' % (key, med, mean, std, nev)

            elif nev < 5:
                single_stat_in = '''
                        <tr>\n
                            <td style="color:#c8c8c8">%s</td>\n
                            <td style="color:#c8c8c8">%s</td>\n
                            <td style="color:#c8c8c8">%s</td>\n
                            <td style="color:#c8c8c8">%s</td>\n
                            <td style="color:#c8c8c8">%s</td>\n
                        </tr>\n''' % (key, med, mean, std, nev)

            elif med > 10 or med < 0.1:
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
            <h4>Amplitude gain ratio test results - statistics</h4>
            <p style="font-size:50%">Grey color: less than 5 results; red: Amplitude ratios larger than 10 or smaller than 0.1.</p>\n'''

        if len(method) == 2:
            string1a = '''
                <p style="font-size:50%">Method: Relative to reference station (green).</p>\n'''
        
        elif method == 'scale_one':
            string1a = '''
                <p style="font-size:50%">Method: Relative to scale of 1.</p>\n'''
        
        elif method == 'median_all_avail':
            string1a = '''
                <p style="font-size:50%">Method: Relative to median of all stations.</p>\n'''

        elif method == 'syn_comp':
            string1a = '''
                <p style="font-size:50%">Method: Relative to synthetic data.</p>\n'''

        string1b = '''
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
        string = string1+string1a+string1b+string2+string3

        return string





def add_single_station_results(comp, gain_result_dir):

    file = '%s/gains_all_events%s.txt' % (gain_result_dir, comp)
    test = glob.glob(file)
    if not test:

        error = '''
        <section>
            <h4>Amplitude gain ratio test results - single events and stations</h4>
            <p style="font-size:100%">File not found.</p>\n
        </section>
        '''
        
        return error

    else:
        with open(file, 'r') as f:
            lines = f.readlines()

            _list = []


            for i_l,l in enumerate(lines[1:]):
                splitted = l.split(',')

                do = True
                if do is True:

                    if i_l == 0:
                        a = '''
                            <thead><tr>
                            '''
                        b = ['<th>%s</th>\n' %x for x in splitted]

                        c = '''
                            </tr></thead>
                            <tbody>
                            '''
                        line_in = a
                        line_b = ''
                        for bb in b:
                            line_b += bb
                        line_in += line_b+ c

                        #print(line_b)
                        #sys.exit()

                    else:
                        b = []
                        for i_x, x in enumerate(splitted):
                            if x != 'nan' and i_x != 0:
                                b.append('<td>%s</td>\n' % round(float(x),3))
                            else:
                                b.append('<td>%s</td>\n' % x)

                        line_in = '''
                        <tr>\n
                        '''
                        for bb in b:
                            line_in += bb
                        line_in += '</tr>\n'

                    _list.append(line_in)
                    line_in=''

                results = ''
                for s in _list:
                    results += s


            string1 = '''
            <section>
                <h4>Amplitude gain test results - single events and stations</h4>
                <div style="overflow:scroll; height:500px; font-size: 40%;">
                <table>
                ''' #% ('40%'.format())
            #print(string1)
            string2 = results
            string3 = '''
                    </tbody>
                </table>
                </div>
            </section>
            '''
            string = string1+string2+string3
            return string


def add_gain_allevents_plot(comp, gain_result_dir):
    gain_plot = os.path.join(gain_result_dir, 'gains_all_events%s_all_gains.png' % comp)

    if not glob.glob(gain_plot):

        logger.info('Error: gain plot not found for %s' % comp)
        error = ''' \t \t \t \t<section data-transition="fade">\n
                    \t \t \t \t \t <h4>All gain test results per station</h4>\n
                    \t \t \t \t \t <p style="font-size:%s"> %s component</p>\n
                    \t \t \t \t \t<p style="font-size:%s">Figure not found.</p>\n
                    \t \t \t \t</section>
                    ''' % ('50%'.format(), comp, '50%'.format())

        return error

    else:
        gain_plot = os.path.join('../results/gains', gain_plot.split('/')[-1])

        include_this = '''
        <section  data-transition="fade">\n
        <h4>All gain test results per station.</h4>\n
        <p style="font-size:%s"> %s component</p>\n
        <div class="map touch-none">
          <img style="width=auto;" src="%s" />\n
        </div>
        </section>
        ''' % ('50%'.format(), comp, gain_plot)

        return include_this