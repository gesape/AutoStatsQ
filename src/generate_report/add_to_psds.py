import glob, os, logging
from pyrocko import guts


logger = logging.getLogger('REPORT-PSDS')


def add_psds_stats(comp, psd_result_dir):
    psds_plots = glob.glob(os.path.join(psd_result_dir, '*__%s.png' % comp))

    if len(psds_plots) < 1:
        logger.info('Error: psd plots not found for %s' % comp)
        error = ''' \t \t \t \t<section data-transition="fade">\n
                    \t \t \t \t \t <h4>Comparison of synthetic and observed PSDs</h4>\n
                    \t \t \t \t \t <h4>%s component</h4>\n
                    \t \t \t \t \t<p style="font-size:%s">Map not found.</p>\n
                    \t \t \t \t</section>
                    ''' % (comp, '50%'.format())

        return error

    else:
        str1 = '''
                <section>
                    <h4>Comparison of synthetic and observed PSDs</h4>
                    <p style="font-size:50%">Red: synthetic; black: observed.</p>\n

                    <div style="overflow:scroll; height:500px; font-size: 30%;">\n'''
        str2 = '''<p style="font-size:%s">%s component</p>\n''' % ('100%'.format(), comp)

        for ot in psds_plots:
            pl = ot.split('/')[-1]
            st = '%s.%s.%s.%s' % (pl.split('_')[0], pl.split('_')[1], pl.split('_')[2], pl.split('_')[3][0])
            str2+='''<p>%s</p>\n''' % st
            str2 += '''<img style="margin:0; margin-bottom:-0.6em; margin-right:0.3em; height: 18em;"
                                src="../results/freq/%s" />\n''' % (pl)
        str3 = '''
                </section>'''

        return str1+str2+str3


def add_psds_flatranges(comp, psd_result_dir):
    psds_plots = glob.glob(os.path.join(psd_result_dir, '*%s_flatrange*.png' % comp))

    if len(psds_plots) < 1:
        logger.info('Error: psd plots flat ranges not found for %s' % comp)
        error = ''' \t \t \t \t<section data-transition="fade">\n
                    \t \t \t \t \t <h4>Flat frequency ranges</h4>\n
                    \t \t \t \t \t <h4>%s component</h4>\n
                    \t \t \t \t \t<p style="font-size:%s">Map not found.</p>\n
                    \t \t \t \t</section>
                    ''' % (comp, '50%'.format())

        return error

    else:

        str1 = '''
                <section>
                    <h4>Flat frequency ranges</h4>
                    <p style="font-size:50%">Median of the ratio of the synthetic and observed PSDs presented above at each frequency step. Dashed lines show the results of the line fits. The red line below indicates a frequency range with a stable ratio between the spectra of observed and synthetic data.</p>\n

                    <div style="overflow:scroll; height:500px; font-size: 50%;">\n'''
        str2 = '''<p style="font-size:%s">%s component</p>\n''' % ('100%'.format(), comp)
        for ot in psds_plots:
            pl = ot.split('/')[-1]
            str2 += '''<img style="margin:0; margin-bottom:-0.6em; margin-right:0.3em; height: 25em;"
                                src="../results/freq/%s" />\n''' % (pl)
        str3 = '''
                </section>'''

        return str1+str2+str3


def include_table_psds_flatranges(psd_result_dir):
    try:
        psd_stats = guts.load(filename='%s/psd_flat_ratio_linefit2.yaml' % psd_result_dir)
    except FileNotFoundError:
        logger.info('Error: PSD flat ranges table not found.')
        error = ''' \t \t \t \t<section data-transition="fade">\n
                    \t \t \t \t \t <h4>PSD test results - flat frequency ranges</h4>\n
                    \t \t \t \t \t<p style="font-size:%s">File not found.</p>\n
                    \t \t \t \t</section>
                    ''' % ('50%'.format())

        return error

    flatranges = psd_stats.FlatFreqRanges
    
    stats_list = []
    for st in sorted(flatranges):

        range_start = round(flatranges[st][0][0],3)
        range_end = round(flatranges[st][0][1],3)
        
        single_stat_in = '''
                <tr>\n
                    <td>%s</td>\n
                    <td>%s</td>\n
                    <td>%s</td>\n
                </tr>\n''' % (st, range_start, range_end)

        stats_list.append(single_stat_in)

    stats = ''
    for s in stats_list:
        stats += s

    string1 = '''
    <section>
        <h4>PSD test results - flat frequency ranges</h4>
        <p style="font-size:50%">Note: Maximum end of the range is defined by the filtering settings and the sampling of the Green's function database which is used to compute the synthetic data. An earlier end (and a later start) points of larger deviations between observed and synthetic data.</p>\n
        <div style="overflow:scroll; height:500px; font-size: 40%;">
        <table>
            <thead><tr>
                <th>Station</th>
                <th>Start [Hz]</th>
                <th>End [Hz]</th>
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

