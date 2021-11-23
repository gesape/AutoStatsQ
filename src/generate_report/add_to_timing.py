import glob, os, logging
from pyrocko import guts
from autostatsq.timing import results_dict

logger = logging.getLogger('REPORT-TIMING')

def add_timing_results(time_result_dir):
    ''' read and add result table to html report'''
    try:
        time_stats = guts.load(filename='%s/timings.yaml' % time_result_dir)
    except FileNotFoundError:
        logger.info('Error: Time test result table not found.')
        error = ''' \t \t \t \t<section data-transition="fade">\n
                    \t \t \t \t \t <h4>Timing test results - statistics</h4>\n
                    \t \t \t \t \t<p style="font-size:%s">File not found.</p>\n
                    \t \t \t \t</section>
                    ''' % ('50%'.format())

        return error

    medians = time_stats.medians
    means = time_stats.means
    stds = time_stats.st_devs
    nevs = time_stats.n_ev

    stats_list = []

    for key in medians.keys():

        med = round(medians[key],1)
        mean = round(means[key],1)
        std = round(stds[key],1)
        nev = nevs[key]

        if nev < 3:
            single_stat_in = '''
                    <tr>\n
                        <td style="color:#c8c8c8">%s</td>\n
                        <td style="color:#c8c8c8">%s</td>\n
                        <td style="color:#c8c8c8">%s</td>\n
                        <td style="color:#c8c8c8">%s</td>\n
                        <td style="color:#c8c8c8">%s</td>\n
                    </tr>\n''' % (key, med, mean, std, nev)

        elif abs(med) > 3:
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
        <h4>Timing test results - statistics</h4>
        <p style="font-size:50%">Grey font color for stations with less than 3 events above cross-correlation threshold, red font color for timing errors larger 3s. Please be aware that this test is rather new and needs more testing.</p>\n
        <div style="overflow:scroll; height:500px; font-size: 40%;">
        <table>
            <thead><tr>
                <th>Station</th>
                <th>Median [s]</th>
                <th>Mean [s]</th>
                <th>Std dev [s]</th>
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


def add_timing_array(time_result_dir):

    timing_plot = os.path.join(time_result_dir, 'timing_arrays.png')

    if not glob.glob(timing_plot):

        logger.info('Error: Timing plot arrays not found.')
        error = ''' \t \t \t \t<section data-transition="fade">\n
                    \t \t \t \t \t <p style="font-size:50%">Array plots showing time shift of maximum cross-correlation value between observed and synthetic data.</p>\n
                    \t \t \t \t \t<p style="font-size:50%">Plot not found.</p>\n
                    \t \t \t \t</section>
                    '''

        return error
    else:
        str1 = '''
                <section>
                    <p style="font-size:50%">Array plots showing time shift of maximum cross-correlation value between observed and synthetic data.</p>
                    <div style="overflow:scroll; height:500px; font-size: 30%;">\n'''

        str2 = '''<img style="margin:0; margin-bottom:-0.6em; margin-right:0.3em; width: %s;"
                                src="../results/timing/timing_arrays.png" />\n''' % ('50%'.format())
        str3 = '''
                </section>'''

        return str1+str2+str3


def add_timing_perstation(time_result_dir):

    timing_plot = os.path.join(time_result_dir, 'timing_errors_allStats.png')

    if not glob.glob(timing_plot):

        logger.info('Error: Timing plot all stations not found.')
        error = ''' \t \t \t \t<section data-transition="fade">\n
                    \t \t \t \t \t <p style="font-size:50%">Time shifts per station; mean and single events.</p>\n
                    \t \t \t \t \t<p style="font-size:50%">Plot not found.</p>\n
                    \t \t \t \t</section>
                    '''

        return error
    else:

        #str1 = '''
        #        <section>
        #            <p style="font-size:50%">Time shifts per station; mean and single events.</p>
        #            <div style="width:200;height:200;overflow:scroll;">\n'''

        #str1 = '''
        #        <section>
        #        <p style="font-size:50%">Time shifts per station; mean and single events.</p>
        #        <div class="box box-three">
        #        '''

        #str1 = '''
        #        <section>
        #        <p style="font-size:50%">Time shifts per station; mean and single events.</p>
        #        <div class="container">
        #        '''
        #str1 = '''
        #        <section>
        #        <p style="font-size:50%">Time shifts per station; mean and single events.</p>
        #        <div id="myDiv" class="img-wrapper">
        #        '''

        #str2 = '''<img class="object-fit-cover" "style="min-width: %s; height: inherit !important;" src="%s" />\n''' % ('100%'.format(), timing_plot)
        #str2 = '''<img src="%s" />\n''' % (timing_plot)

        str1 = '''
        <section>
        <div class="map touch-none">
          <img style="height: 5000px, width=auto;" src="../results/timing/timing_errors_allStats.png" />\n
        '''

        str3 = '''
                </div>
                </section>'''

        return str1+str3