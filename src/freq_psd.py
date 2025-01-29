import numpy as num
import math
import os
import logging

import matplotlib.pyplot as plt
from pyrocko import util, pile, trace
from pyrocko.guts import Object, Dict, String, Float, List, Tuple


deg = 1


def flat_by_neighbor_comp(mean_rat, cha, st, f_syn_keep,
                          dir_f, n, fac_norm, f_ign,
                          plot_psd_neighbcomp):
    """
    input: mean psd ratio
    ...currently not used!...
    ...testing different methods to get flat f range...

    calc difference of each point to next point (or to prev and next?)
    plot that!
    --> then get flat range (look at data, maybe obvious, maybe line fit?)
    """

    diff_next = [abs(mean_rat[i+1] - mean_rat[i])
                 for i in range(len(mean_rat[0:-1]))]

    mean_diff_next = num.median(diff_next)
    diff_next_rel = [dn / mean_diff_next for dn in diff_next]

    if plot_psd_neighbcomp is True:
        fig, ax = plt.subplots(nrows=8, ncols=1, figsize=(12, 10), sharex=True)
        ax[0].plot(f_syn_keep, mean_rat, '.')
        # plot diff to next
        ax[1].plot(f_syn_keep[0:-1], diff_next, '.')
        ax[1].set_title('abs. dist to next neighbour')
        # plot diff to next / median(mean) of diffs
        ax[2].plot(f_syn_keep[0:-1], diff_next_rel, '.')
        ax[2].set_title('abs. dist/ mean of diffs')
        # plot range for which diff/mean is smaller than threshold
        for i_df, df in enumerate(diff_next_rel):
            if df < 0.5:
                ax[3].plot(f_syn_keep[i_df], 1, 'b.')
        ax[3].set_title('threshold 0.1')
        ax[3].set_yticks([])
        # plot sum of abs diffs for 5 f points
    sum_diff = []

    for i_df, df in enumerate(diff_next[:-5]):
        sum_diff.append(sum(diff_next[i_df:i_df+5]))

    if plot_psd_neighbcomp is True:
        ax[4].plot(f_syn_keep[:-6], sum_diff, '.')
        ax[4].set_title('summed diffs 5 pts')
        # plot polyfit through summed diffs

    n = 15
    flat_by_sumdiff = []
    flat_by_sumdiff_y = []
    go_in = True
    norm_reg = 100
    lsq_sqthresh = 30
    dev = int(num.floor(n/2))

    for i in range(int(num.floor(n/2)),
                   int(len(sum_diff)-num.floor(n/2))):

        startp = int(i-dev)
        stopp = int(i+dev)

        x = f_syn_keep[startp:stopp]
        y = sum_diff[startp:stopp]

        reg, res_lsq, _, _, _ = num.polyfit(x=x, y=y, deg=deg, full=True)

        p = num.poly1d(reg)
        xp = num.linspace(f_syn_keep[startp], f_syn_keep[stopp], 100)

        if plot_psd_neighbcomp is True:
            ax[5].plot(xp, p(xp))
            ax[6].plot(f_syn_keep[i-dev], num.sqrt(res_lsq/n), 'b.')

        if abs(reg[0])/norm_reg < 1 and\
          num.sqrt(res_lsq/n) < lsq_sqthresh:

            if go_in is True:
                f_start = f_syn_keep[i-dev]
                f_stop = f_syn_keep[i+dev]
                go_in = False

            if f_syn_keep[i-dev] < f_stop or\
              f_syn_keep[i-dev] < f_stop + f_ign:
                f_stop = f_syn_keep[i+dev]

            if f_syn_keep[i-dev] > f_stop:
                flat_by_sumdiff.append((f_start, f_stop))
                flat_by_sumdiff_y.append(num.nanmean(y))
                f_start = f_syn_keep[i-dev]
                f_stop = f_syn_keep[i+dev]

    flat_by_sumdiff.append((f_start, f_stop))
    flat_by_sumdiff_y.append(num.nanmean(y))

    if plot_psd_neighbcomp is True:
        for frange in flat_by_sumdiff:
            y = [1, 1]
            ax[7].plot(frange, y, '-b')
        ax[5].set_title('line fit through summed diffs, only if < threshold')
        ax[6].set_title('sqrt(lsq/n_point)')
        ax[7].set_title('*flat* acc. to slope and lsq threshold')
        ax[7].set_yticks([])

        ax[0].set_title(str(st.station)+cha)
        plt.tight_layout()
        fig.savefig(os.path.join(dir_f, '%s_%s_%s_flattests.png' % (st.network, st.station, cha)))
        plt.close(fig)

    return flat_by_sumdiff, flat_by_sumdiff_y


def calc_reg_prepplot(npoints, i, f, r, deg=1, plot=False):
    """
    Calculate regression of degree deg over n points and
    return data vectors for plotting.

    :param npoints: number of points to use for fitting line
    :paran i: current index in for-loop
    :param deg: polynomial degree, default is 1
    :param f: frequencies (1-D np array or list)
    :param r: amplitude data (1-D np array or list)
    """
    n = npoints
    startp = int(i-n/2)
    stopp = int(i+n/2)

    x = f[startp:stopp]
    y = r[startp:stopp]

    reg, res_lsq, _, _, _ = num.polyfit(x=x, y=y, deg=deg, full=True)

    if plot is True:
        p = num.poly1d(reg)
        xp = num.linspace(f[startp], f[stopp], 100)

        return reg, res_lsq, p, xp
    else:
        return reg, res_lsq


def get_flat_freq_ranges(r, f, n, fac_norm, f_ign, flat_areas):
    """
    get flat frequency ranges from indexes provided in flat_areas list
    called by get_flat_freqs

    for single nslc spec ratio
    """
    _flat_f_ranges = []
    flat_ranges_all = []
    y_flat_f_ranges = []

    for i, no in enumerate(flat_areas[0:-1]):
        flat_ranges_all.append((f[no-int(n/2)], f[no+int(n/2)]))
        y_flat_f_ranges.append(r[no-int(n/2):no+int(n/2)])

    flat_ranges_sorted = []
    y_flat_f_ranges_sorted = []

    for i_tup, tup in enumerate(flat_ranges_all):
        if i_tup == 0:
            current_start = tup[0]
            current_stop = tup[1]
            i_tup_y_start = 0

        else:
            if tup[0] < current_stop + f_ign:
                current_stop = tup[1]
            else:
                flat_ranges_sorted.append((current_start, current_stop))
                current_start = tup[0]
                current_stop = tup[1]
                y_flat_f_ranges_sorted.append(
                    num.median(
                        y_flat_f_ranges[i_tup_y_start:i_tup]))
                i_tup_y_start = i_tup

        if i_tup == len(flat_ranges_all)-1:
            flat_ranges_sorted.append((current_start, current_stop))
            y_flat_f_ranges_sorted.append(
                num.median(
                    y_flat_f_ranges[i_tup_y_start:i_tup]))

    _flat_f_ranges = flat_ranges_sorted

    return _flat_f_ranges, y_flat_f_ranges_sorted


def get_flat_freqs(n, freqs, rat_a, fac_norm, f_ign,
                   deg=1, plot_flat_range=False):
    """
    get frequency ranges at which the ratio of syn and obs spect are flat
    for single nslc spec ratio

    :param n: number of points to use for fitting line
    :param freqs: frequencies (1-D np array or list)
    :param rat_a: amplitude ratio data (1-D np array or list)
    :param fac_norm: normalization factor for slope
    :param deg: polynomial degree, default is 1

    """

    r = rat_a
    f = freqs
    # df = f[1]-f[0]

    # plot data
    if plot_flat_range is True:
        fig4, ax4 = plt.subplots(nrows=2, ncols=1, figsize=(8, 6), sharex=True,
                                 gridspec_kw={'height_ratios': [6, 1]})

        ax4[0].plot(f, r, 'k.', markersize=12)
        ax4[0].set_ylabel('median PSD ratio (syn/obs)')

    rmin = num.min(r)
    rmax = num.max(r)

    flat_areas = []
    lsq_sqthresh = 10
    dev = int(num.floor(n/2))

    for i in range(dev, int(len(r)-dev)):
        if plot_flat_range is True:
            reg, res_lsq, p, xp = calc_reg_prepplot(n, i, f, r, deg,
                                                    plot=plot_flat_range)
            ax4[0].plot(xp, p(xp), '--')
            ax4[0].set_ylim(0, 1*fac_norm)
        else:
            reg, res_lsq = calc_reg_prepplot(n, i, f, r, deg,
                                             plot=plot_flat_range)

        if abs(reg[0])/fac_norm < 1 and num.sqrt(res_lsq/n) < lsq_sqthresh:
            flat_areas.append(i)

    _flat_f_ranges, y_flat_f_ranges_sorted =\
        get_flat_freq_ranges(r, f, n, fac_norm, f_ign, flat_areas)

    if plot_flat_range:
        return _flat_f_ranges, y_flat_f_ranges_sorted, rmin, rmax, fig4, ax4
    else:
        return _flat_f_ranges, y_flat_f_ranges_sorted, rmin, rmax


class dict_stats(Object):
    """
    Dict for all stations + their flat freq ranges
    """
    FlatFreqRanges = Dict.T(String.T(), List.T(Tuple.T(2, Float.T())))
    MeanMedianR_FlatRanges = Dict.T(String.T(), List.T(Float.T()))


def dump_flat_ranges(flat_f_ranges_stlist, freq_rat_list_y,
                     nslc_list, dir_f, fname_ext, only_first=True):
    """
    writing flat areas of all stations to file

    """
    logs = logging.getLogger('dump_flat_ranges')
    logs.info('start writing to file')

    f_r = []
    f_y = []
    if only_first is True:
        f_r = [[f[0]] if len(f) > 1 else f for f in flat_f_ranges_stlist]
        f_y = [[f[0]] if len(f) > 1 else f for f in freq_rat_list_y]

    else:
        f_r = flat_f_ranges_stlist
        f_y = freq_rat_list_y

    fr_dict = dict(zip(map(lambda x: x, nslc_list), f_r))
    y_dict = dict(zip(map(lambda x: x, nslc_list), f_y))

    flfr = dict_stats(FlatFreqRanges=fr_dict,
                      MeanMedianR_FlatRanges=y_dict)
    flfr.regularize()
    flfr.validate()

    flfr.dump(filename=os.path.join(dir_f, 'psd_flat_ratio_%s.yaml' % fname_ext))


def const_psd_rat(mean_rat, cha, st, l, f_syn_keep,
                  n, fac_norm, f_ign,
                  dir_f=False, plot_flat_range=False):
    """
    Find frequency range of roughly constant relative gain factor
    (optimal is no needed gain factor) in mean psd ratio of syn and
    obs data.

    :param mean_rat: Mean of PSD ratios for current nslc over all events
    """

    net = st.network
    stat = st.station
    # num.save('mn_rats_%s_%s_%s.np' % (net, stat, cha), mean_rat)
    # num.save('freqs', f_syn_keep)

    # data
    freqs = f_syn_keep
    rat_a = mean_rat
    assert len(freqs) == len(rat_a)

    if plot_flat_range:
        _flat_f_ranges, y_flat_f_ranges_sorted, rmin, rmax, fig4, ax4 =\
        get_flat_freqs(n, freqs, rat_a, fac_norm, f_ign, deg=deg,
                       plot_flat_range=plot_flat_range)

    else:
        _flat_f_ranges, y_flat_f_ranges_sorted, rmin, rmax = get_flat_freqs(
            n, freqs, rat_a, fac_norm, f_ign, deg=1)

    flat_f_ranges = [i for i in _flat_f_ranges if i[0] != i[1]]
    y_flat_f_ranges = [y for y, i in zip(
                       y_flat_f_ranges_sorted, _flat_f_ranges)
                       if i[0] != i[1]]

    if plot_flat_range is True:
        ax4[0].set_title('%s %s %s' % (net, stat, cha))

        for fl in flat_f_ranges:
            l_fr = num.linspace(fl[0], fl[1], 100)
            l_y = [1 for i in range(len(l_fr))]
            ax4[1].plot(l_fr, l_y, 'r', markersize=2)
            ax4[1].yaxis.set_ticks([])
            ax4[1].set_xlabel('Frequency [Hz]')

        fig4.savefig(os.path.join(dir_f, '%s_%s_%s_%s_flatrange%s_%s.png' % (net, stat, l, cha, n, fac_norm)))

        fig4.tight_layout()
        plt.close(fig4)

    return flat_f_ranges, y_flat_f_ranges


def plot_m_ratio(median_rat, f_syn_Z_nonan, nsl, l, dir_f, cha):
    """
    Plotting of mean ratio of syn and obs PSD over all events
    (for current channel).

    """
    fig3, ax3 = plt.subplots(nrows=1, ncols=1, figsize=(6, 5))

    if not median_rat == []:
        ax3.plot(f_syn_Z_nonan[1:], median_rat[1:], 'k.')
        fmax = max(f_syn_Z_nonan)
        fmin = 0.01#0.001
        ax3.set_xlim(fmin, fmax)
        #ax3.set_xscale('log')

        ax3.set_xlabel('Frequency [Hz]')
        ax3.set_ylabel('median PSD ratio (syn/obs)')

        # ax3.set_yscale('log')
        fig3.tight_layout()
        fig3.savefig(os.path.join(dir_f, '%s_%s_%s_%s_mratio.png' % (nsl[0], nsl[1], l, cha)))
        plt.close(fig3)


def calc_mean_ratio(ratio_npar):   # hier dran arbeiten!
    '''
    Calculate mean of ratios for one station over all events.
    '''
    mean_rat = num.nanmean(ratio_npar, axis=0)
    return(mean_rat)


def calc_median_ratio(ratio_npar):   # hier dran arbeiten!
    '''
    Calculate mean of ratios for one station over all events.
    '''
    median_rat = num.nanmedian(ratio_npar, axis=0)
    return(median_rat)


def plot_psdratio_from_dict(ratpsd_by_event, st, l, cha, catalog, dir_f):
    n_ev = len(catalog)
    nrows = math.ceil(n_ev / 5)
    ncols = 5
    fmin = 0.01#0.001
    fmax = None

    r_max = 0
    r_min = 0
    for k, i in ratpsd_by_event.items():
        if not i[1] == []:
            r_maxnow = num.max(i[1])
            r_minnow = num.min(i[1])
            if r_maxnow > r_max:
                r_max = r_maxnow
            if r_minnow < r_min:
                r_min = r_minnow

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(13, nrows*3))

    for i_ev, ev in enumerate(catalog):

        if nrows > 1:
            i_x = int(i_ev/ncols)
            i_y = int(i_ev % ncols)

            if i_x == nrows-1:
                ax[i_x, i_y].set_xlabel('Frequency [Hz]')

            if i_y == 0:
                ax[i_x, i_y].set_ylabel('PSD ratio (syn/obs)')

            ev_time_str = util.time_to_str(ev.time)[0:10]
            ax[i_x, i_y].set_title(ev_time_str)
            #ax[i_x, i_y].set_xscale('log')
            #ax[i_x, i_y].set_yscale('log')

            try:
                f_syn, rat = ratpsd_by_event[ev_time_str]
                ax[i_x, i_y].plot(f_syn[1:], rat[1:], 'k.')
                if i_ev == 0 or not fmax:
                    fmax = max(f_syn)
            except ValueError:
                logging.error('ValueError obs')
            except KeyError:
                logging.error('no data key obs')

            ax[i_x, i_y].set_xlim(fmin, fmax)
            ax[i_x, i_y].set_ylim(0, 10)

        if nrows*ncols > n_ev:
            dif = nrows*ncols - n_ev
            for i in range(dif):
                ax[i_x, i_y+i+1].set_xlabel('Frequency [Hz]')

        if nrows == 1:
            i_x = int(i_ev)
            ax[i_x].set_xlabel('Frequency [Hz]')

            ax[i_x].set_ylabel('PSD ratio (syn/obs)')

            ev_time_str = util.time_to_str(ev.time)[0:10]
            ax[i_x].set_title(ev_time_str)
            #ax[i_x, i_y].set_xscale('log')
            #ax[i_x, i_y].set_yscale('log')

            try:
                f_syn, rat = ratpsd_by_event[ev_time_str]
                ax[i_x].plot(f_syn[1:], rat[1:], 'k.')
                if i_ev == 0 or not fmax:
                    fmax = max(f_syn)
            except ValueError:
                logging.error('ValueError obs')
            except KeyError:
                logging.error('no data key obs')

            ax[i_x].set_xlim(fmin, fmax)
            ax[i_x].set_ylim(0, 10)
    try:
        plt.tight_layout()
    except Exception:
        pass
    try:
        #print('here')
        fig.savefig(os.path.join(dir_f, '%s_%s_%s_%s_ratio.png' % (st.network, st.station, l, cha)))
    except Exception:
        pass
    # plt.show()
    plt.close(fig)


def plot_psd_from_dict(obspsd_by_event, synpsd_by_event,
                       st, l, cha, catalog, dir_f):

    n_ev = len(catalog)
    nrows = math.ceil(n_ev / 5)
    ncols = 5
    fmin = 0.01#0.001
    fmax = None

    a_max = 0
    a_min = 0
    for k, i in obspsd_by_event.items():
        if not i[1] == []:
            a_maxnow = num.max(i[1])
            a_minnow = num.min(i[1])
            if a_maxnow > a_max:
                a_max = a_maxnow
            if a_minnow < a_min:
                a_min = a_minnow

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(13, nrows*3))

    for i_ev, ev in enumerate(catalog):
        if nrows > 1:
            i_x = int(i_ev/ncols)
            i_y = int(i_ev % ncols)

            if i_x == nrows-1:
                ax[i_x, i_y].set_xlabel('Frequency [Hz]')

            if i_y == 0:
                ax[i_x, i_y].set_ylabel('PSD')

            ev_time_str = util.time_to_str(ev.time)[0:10]

            ax[i_x, i_y].set_title(ev_time_str)
            ax[i_x, i_y].set_xscale('log')
            ax[i_x, i_y].set_yscale('log')

            try:
                f_obs_Z, a_obs_Z = obspsd_by_event[ev_time_str]
                ax[i_x, i_y].plot(f_obs_Z[1:], a_obs_Z[1:], 'k')
                if i_ev == 0 or not fmax:
                    fmax = max(f_obs_Z)
            except ValueError:
                logging.error('ValueError obs')
            except KeyError:
                logging.error('no data key obs: %s' % ev_time_str)
            try:
                f_syn_Z, a_syn_Z = synpsd_by_event[ev_time_str]
                ax[i_x, i_y].plot(f_syn_Z[1:], a_syn_Z[1:], 'r')
                if i_ev == 0 or not fmax:
                    fmax = max(f_syn_Z)
            except ValueError:
                logging.error('ValueError syn')
            except KeyError:
                logging.error('no data key syn: %s' % ev_time_str)

            ax[i_x, i_y].set_xlim(fmin, fmax)
            ax[i_x, i_y].set_ylim(a_min, a_max)

        if nrows*ncols > n_ev:
            dif = nrows*ncols - n_ev
            for i in range(dif):
                ax[i_x, i_y+i+1].set_xlabel('Frequency [Hz]')

        if nrows == 1:
            i_x = i_ev

            ax[i_x].set_xlabel('Frequency [Hz]')
            ax[i_x].set_ylabel('PSD')

            ev_time_str = util.time_to_str(ev.time)[0:10]

            ax[i_x].set_title(ev_time_str)
            ax[i_x].set_xscale('log')
            ax[i_x].set_yscale('log')

            try:
                f_obs_Z, a_obs_Z = obspsd_by_event[ev_time_str]
                ax[i_x].plot(f_obs_Z[1:], a_obs_Z[1:], 'k')
                if i_ev == 0 or not fmax:
                    fmax = max(f_obs_Z)
            except ValueError:
                logging.error('ValueError obs')
            except KeyError:
                logging.error('no data key obs: %s' % ev_time_str)
            try:
                f_syn_Z, a_syn_Z = synpsd_by_event[ev_time_str]
                ax[i_x].plot(f_syn_Z[1:], a_syn_Z[1:], 'r')
                if i_ev == 0 or not fmax:
                    fmax = max(f_syn_Z)
            except ValueError:
                logging.error('ValueError syn')
            except KeyError:
                logging.error('no data key syn: %s' % ev_time_str)

            ax[i_x].set_xlim(fmin, fmax)
            ax[i_x].set_ylim(a_min, a_max)

    plt.tight_layout()
    fig.savefig(os.path.join(dir_f, '%s_%s_%s_%s.png' % (st.network, st.station, l, cha)))
    plt.close(fig)


def get_a_f(traces, cha):
    """
    Calculation of actual PSD (freqs and ampls).

    :returns: 1-D numpy arrays f and a with frequencies and amplitudes.
    """

    a_list = []

    for tr in traces:
        # print(tr.channel)
        # if tr.channel not in ('Z', 'R','T'):
        #    continue
        win = num.ones(tr.data_len())
        # tr.snuffle()

        # demean:
        tr.ydata -= tr.ydata.mean()
        tr.ydata *= win

        # calc spectrum:
        f, a = tr.spectrum(pad_to_pow2=True)
        a = num.abs(a)**2   # square ampls
        a *= tr.deltat * 2. / num.sum(win**2)  # normalizing
        a[0] /= 2.  # erster wert durch zwei geteilt.
        a[a.size//2] /= 2.  # mittlerer wert durch zwei geteilt.
        a_list.append(a)

    # print(len(a_list))
    # print([len(a) for a in a_list])

    try:
        a = num.vstack(a_list)
        a = num.median(a, axis=0)

        return f, a

    except ValueError:
        # print('ValueError %s' %cha)
        f = []
        a = []
        return f, a


def calc_plot_psds(catalog, data_pile, syn_data_pile,
                   cha, l, dir_f, arrT_array, arrT_R_array, 
                   nsl, i_st, nst,
                   tinc, tpad, dt_s, dt_e,
                   n, fac_norm, f_ign,
                   plot_psds, plot_ratio_extra,
                   plot_m_rat, plot_flat_ranges): # ,
                   # plot_neighb_ranges):
    """
    Next level function, called by ```prep_psd_fct``` to procede with current
    synthetic and observed datapiles.

    Selects for each channel and event traces from the piles and passes them
    on to function ```get_a_f``` that calculates the PSD.

    For parameters see ```prep_psd_fct```.

    """

    n_ev = len(catalog)
    ratio_npar = num.empty((1,))
    f_syn_keep = 0

    #if plot_psds is True:
    obspsd_by_event = {}
    synpsd_by_event = {}

    #if plot_ratio_extra is True:
    ratpsd_by_event = {}

    cnt = 0

    for i_ev, ev in enumerate(catalog):
        print('Station: %5d/%s, Event: %5d/%s' % (i_st, nst, i_ev,n_ev), end='\r')
        ev_time_str = util.time_to_str(ev.time)[0:10]
        # print(ev_time_str)
        # abs. arrival times: arrT_array,
        # shape (len(subset_catalog), len(ns))

        start_twd = arrT_array[i_ev, i_st] - dt_s
        end_twd = arrT_R_array[i_ev, i_st] - dt_e

        #print(util.time_to_str(start_twd))
        #print(util.time_to_str(end_twd))


        trs_obs = data_pile.all(
            tmin=start_twd,
            tmax=end_twd,
            tinc=tinc,
            tpad=tpad,
            trace_selector=lambda tr: tr.nslc_id[:1] == nsl[:1] and
                                      tr.nslc_id[3] == cha
                                      and tr.nslc_id[2] == l,
            want_incomplete=True)

        trs_syn = syn_data_pile.all(
            tmin=start_twd,
            tmax=end_twd,
            tinc=tinc,
            tpad=tpad,
            trace_selector=lambda tr: tr.nslc_id[:1] == nsl[:1] and
                                      tr.nslc_id[3] == cha,
            want_incomplete=True)


        if trs_syn and trs_obs:
            #print(len(trs_obs[0].ydata), len(trs_syn[0].ydata))
            #trace.snuffle([trs_syn[0], trs_obs[0]])
            #trace.snuffle(trs_obs + trs_syn)
            f_obs_Z, a_obs_Z = get_a_f(trs_obs, cha)
            f_syn_Z, a_syn_Z = get_a_f(trs_syn, cha)

            if plot_psds is True:
                obspsd_by_event['%s' % (ev_time_str)] = (f_obs_Z, a_obs_Z)
                synpsd_by_event['%s' % (ev_time_str)] = (f_syn_Z, a_syn_Z)

            if num.array_equal(f_syn_Z, f_obs_Z) and len(a_syn_Z) > 0:
                # diff_a = abs(a_syn_Z - a_obs_Z)

                ratio_a = [a1/a2 if not a2 == 0 else 0
                           for a1, a2 in zip(a_syn_Z, a_obs_Z)]

                # print([(f, a1, a2, a1/a2, d, '\n') for f,a1,a2,d
                #        in zip(f_obs_Z, a_obs_Z, a_syn_Z, ratio_a)])
                if plot_ratio_extra is True:
                    ratpsd_by_event['%s' % (ev_time_str)] = (f_syn_Z, ratio_a)

                if cnt == 0:
                    ratio_npar = num.empty([n_ev, len(f_syn_Z)])
                    ratio_npar[:] = num.nan
                try:
                    ratio_npar[i_ev, :] = ratio_a
                    cnt += 1
                    f_syn_keep = f_syn_Z
                except ValueError:
                    logging.error('not same length')
    #if plot_psds is True and plot_ratio_extra is True and f_syn_keep != 0:
    return obspsd_by_event, synpsd_by_event, ratpsd_by_event,\
               f_syn_keep, ratio_npar

    # elif plot_psds is True and plot_ratio_extra is not True and f_syn_keep != 0:
    #     return obspsd_by_event, synpsd_by_event, f_syn_keep, ratio_npar

    # elif plot_psds is not True and plot_ratio_extra is True and f_syn_keep != 0:
    #     return ratpsd_by_event, f_syn_keep, ratio_npar

    # else: #plot_psds is not True and plot_ratio_extra is not True:
    #     return f_syn_keep, ratio_npar


def prep_psd_fct(i_st, st, nst, l, subset_catalog, dir_f, arrT_array, arrT_R_array,
                 datapath,
                 syndatapath, tinc, tpad, dt_s, dt_e,
                 n, fac_norm, f_ign,
                 plot_psds=False, plot_ratio_extra=False,
                 plot_m_rat=False, plot_flat_ranges=False): # ,
                 # plot_neighb_ranges=False):
    """
    Preparing the PSD calculations and plotting, e.g. making data-piles
    and calling next function if both, synthetic and recorded data is
    in piles.

    :param i_st: number of station - must fit arrival time table!
    :param st: pyrocko station object (current station)
    :param subset_catalog: Catalog to be used (list of pyrocko events)
    :param dir_f: directory to store PSD results at
    :param arrT_array: Arrival times
    :param datapath: Path to restituted (,rotated and downsampled) data
    :param syndatapath: Path to sythetic data
    :param plot_ratio_extra: Should ratios of PSDs be plotted to extra figure?
                             They are displayed within the PSD plot anyway,
                             default is False.

    :returns freq_rat_list: List containing freuqncy ranges at which the
                            ratio between synthetic and observed psd is
                            constant
    :returns nslc_list: nslc list for current station
    """

    st_name = st.station
    st_data_pile = pile.make_pile(datapath, regex='%s_%s_' % (st.network, st_name),
                                  show_progress=False)
    st_syn_data_pile = pile.make_pile(syndatapath, regex='%s_%s_' % (st.network, st_name),
                                      show_progress=False)
    freq_rat_list = []
    freq_rat_list_y = []
    nslc_list = []
    # flat_by_next = []
    # flat_by_next_y = []

    if st_data_pile.tmin is not None and st_data_pile.tmax is not None and\
        st_syn_data_pile.tmin is not None\
        and st_syn_data_pile.tmax is not None:

        nsl = st.nsl()

        for cha in ['Z', 'R', 'T']:

            outs = calc_plot_psds(
                subset_catalog,
                st_data_pile,
                st_syn_data_pile,
                cha, l, dir_f,
                arrT_array, arrT_R_array,
                nsl, i_st, nst,
                tinc, tpad,
                dt_s, dt_e,
                n, fac_norm, f_ign,
                plot_psds=plot_psds,
                plot_ratio_extra=plot_ratio_extra,
                plot_m_rat=plot_m_rat,
                plot_flat_ranges=plot_flat_ranges)  # ,
                # plot_neighb_ranges=plot_neighb_ranges)

            ratio_npar = outs[-1]
            f_syn_keep = outs[-2]

            if plot_psds is True and outs[0] and outs[1]:
                plot_psd_from_dict(outs[0], outs[1], st, l, cha,
                                   subset_catalog, dir_f)

                if plot_ratio_extra is True and outs[2]:
                    plot_psdratio_from_dict(outs[2], st, l, cha,
                                            subset_catalog, dir_f)

            elif plot_psds is not True and plot_ratio_extra is True and outs[2]:
                plot_psdratio_from_dict(outs[2], st, l, cha,
                                        subset_catalog, dir_f)

            if ratio_npar.shape != (1,):
                m_rat = calc_median_ratio(ratio_npar)

                if plot_m_rat is True:
                    plot_m_ratio(m_rat, f_syn_keep, nsl, l, dir_f, cha)

                f, y = const_psd_rat(m_rat, cha, st, l, f_syn_keep,
                                     n, fac_norm, f_ign,
                                     plot_flat_range=plot_flat_ranges,
                                     dir_f=dir_f)

                freq_rat_list.append(f)
                freq_rat_list_y.append(y)
                nslc_list.append('%s.%s.%s.%s' % (nsl[0], nsl[1], l, cha))
            '''
            f2, y2 = flat_by_neighbor_comp(
                m_rat, cha, st, f_syn_keep, dir_f,
                n, fac_norm, f_ign,
                plot_psd_neighbcomp=plot_neighb_ranges)
            flat_by_next.append(f2)
            flat_by_next_y.append(y2)
            '''

    return freq_rat_list, freq_rat_list_y, nslc_list  # ,\
            # flat_by_next, flat_by_next_y
