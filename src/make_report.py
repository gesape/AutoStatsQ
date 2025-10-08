import logging, os
from shutil import copytree, rmtree

import autostatsq.generate_report.add_to_orient as ori
import autostatsq.generate_report.add_to_cat as cat
import autostatsq.generate_report.add_to_gains as gains
import autostatsq.generate_report.add_to_psds as psds
import autostatsq.generate_report.add_to_timing as time

from autostatsq.config import AutoStatsQConfig


logs = logging.getLogger('REPORT')
logger = logs

file_path = os.path.abspath(os.path.dirname(__file__))

def gen_report(config):
    logger.info('yes')
    tpl_file = os.path.join(file_path, 'generate_report/index_template.html')

    configfile = config  # '../../../../AutoStatsQ_settings.config'

    gensettings, catalogconf, arrTconf, metaDataconf, RestDownconf,\
            synthsconf, gainconf, psdsconf, orientconf, timingconf, tc, maps =\
            AutoStatsQConfig.load(filename=configfile).Settings

    try:
        os.makedirs(os.path.join(gensettings.work_dir, 'result_report'))
    except:
        rmtree(os.path.join(gensettings.work_dir, 'result_report'))
        os.makedirs(os.path.join(gensettings.work_dir, 'result_report'))

    html_directories = ['figures', 'reveal', 'theme']
    for h in html_directories:
        copytree(src=os.path.join(file_path, 'generate_report', h),
                        dst=os.path.join(os.path.join(gensettings.work_dir, 'result_report'), h))

    outfile = os.path.join(gensettings.work_dir, 'result_report', 'index_report.html')

    result_dir = os.path.join(gensettings.work_dir, 'results')

    tmin = catalogconf.tmin_str
    tmax = catalogconf.tmax_str
    Mmin = catalogconf.min_mag

    with open(tpl_file, 'r') as f:
        lines = f.readlines()
        extended_lines = []

        for line in lines:
            if "INCLUDE_HERE_CAT" in line:
                logger.info('Including catalog in html report.')
                cat_result_dir = os.path.join(result_dir, 'catalog')
                extended_lines.append(cat.include_maps_subsets(tmin, tmax, cat_result_dir))
                extended_lines.append(cat.include_event_tabels(tmin, tmax, Mmin, cat_result_dir))
                extended_lines.append(cat.include_statistic_plots(cat_result_dir))

            elif "INCLUDE_HERE_ORIENT" in line:
                logger.info('Including orientation test results in html report.')
                orient_result_dir = os.path.join(result_dir, 'orient')
                extended_lines.append(ori.include_map(orient_result_dir,
                    orientconf.ccmin))
                extended_lines.append(ori.include_table_correction_angles(orient_result_dir,
                    orientconf.ccmin))
                extended_lines.append(ori.add_single_station_results(orient_result_dir,
                    orientconf.ccmin))
                extended_lines.append(ori.add_all_stats_figures_baz(orient_result_dir,
                    orientconf.ccmin))
                extended_lines.append(ori.add_all_stats_figures_timing(orient_result_dir,
                    orientconf.ccmin))
                extended_lines.append(ori.add_all_stats_figures_ccs(orient_result_dir,
                    orientconf.ccmin))
            
            elif "INCLUDE_HERE_GAINS" in line:
                logger.info('Including gain test results in html report.')
                gain_result_dir = os.path.join(result_dir, 'gains')
                
                for comp in ['Z', 'R', 'T']:
                    extended_lines.append(gains.add_gain_map(comp, gain_result_dir))                    
                    extended_lines.append(gains.add_gain_median_mean(comp, gainconf.gain_factor_method, gain_result_dir))
                    extended_lines.append(gains.add_gain_allevents_plot(comp, gain_result_dir))
                    extended_lines.append(gains.add_single_station_results(comp, gain_result_dir))
            
            elif "INCLUDE_HERE_PSDS" in line:
                logger.info('Including PSD test results in html report.')
                psd_result_dir = os.path.join(result_dir, 'freq')
                
                for comp in ['Z', 'R', 'T']:
                    extended_lines.append(psds.add_psds_stats(comp, psd_result_dir))
                    extended_lines.append(psds.add_psds_flatranges(comp, psd_result_dir))

                extended_lines.append(psds.include_table_psds_flatranges(psd_result_dir))

            elif "INCLUDE_HERE_TIMING" in line:
                logger.info('Including timing test results in html report.')
                time_result_dir = os.path.join(result_dir, 'timing')
                extended_lines.append(time.add_timing_results(time_result_dir))
                extended_lines.append(time.add_timing_array(time_result_dir))
                extended_lines.append(time.add_timing_perstation(time_result_dir))

            else:
                extended_lines.append(line)

        with open(outfile, 'w') as o:
            for el in extended_lines:
                o.write(el)

    logger.info('Generated html report and saved in work directory.' +
                ' Use firefox to open index_report.html')
