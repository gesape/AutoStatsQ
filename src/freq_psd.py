import numpy as num
import math

import matplotlib.pyplot as plt
from pyrocko import util, pile
from pyrocko.guts import Object, Dict, String, Float, List, Tuple


def calc_reg_prepplot(npoints, i, f, r, deg=1, plot=False):
    '''
    Calculate regression of degree deg over n points and
    return data vectors for plotting.

    :param npoints: number of points to use for fitting line
    :paran i: current index in for-loop
    :param deg: polynomial degree, default is 1
    :param f: frequencies (1-D np array or list)
    :param r: amplitude data (1-D np array or list)
    '''
    n = npoints
    startp = int(i-n/2)
    stopp = int(i+n/2)

    x = f[startp:stopp]
    y = r[startp:stopp]

    reg = num.polyfit(x=x, y=y, deg=deg)

    if plot is True:
        p = num.poly1d(reg)
        xp = num.linspace(f[startp], f[stopp], 100)

        return reg, p, xp
    else:
        return reg


def get_flat_freq_ranges(r, f, n, flat_areas):
    '''
    get flat frequency ranges from indexes provided in flat_areas list
    called by get_flat_freqs

    for single nslc spec ratio
    '''
    _flat_f_ranges = []
    flat_ranges_all = []

    for i, no in enumerate(flat_areas[0:-1]):
        flat_ranges_all.append((f[no-int(n/2)], f[no+int(n/2)]))

    flat_ranges_sorted = []

    for i_tup, tup in enumerate(flat_ranges_all):
        if i_tup == 0:
            current_start = tup[0]
            current_stop = tup[1]
        else:
            if tup[0] < current_stop:
                current_stop = tup[1]
            else:
                flat_ranges_sorted.append((current_start, current_stop))
                current_start = tup[0]
                current_stop = tup[1]
        if i_tup == len(flat_ranges_all)-1:
            flat_ranges_sorted.append((current_start, current_stop))

    _flat_f_ranges = flat_ranges_sorted

    return _flat_f_ranges


def get_flat_freqs(n, freqs, rat_a, deg=1, plot_flat_range=False):
    '''
    get frequency ranges at which the ratio of syn and obs spect are flat
    for single nslc spec ratio

    :param n: number of points to use for fitting line
    :param freqs: frequencies (1-D np array or list)
    :param rat_a: amplitude ratio data (1-D np array or list)
    :param fac_norm: normalization factor for slope
    :param deg: polynomial degree, default is 1

    '''

    r = rat_a
    f = freqs
    df = f[1]-f[0]

    # plot data
    if plot_flat_range is True:
        fig4, ax4 = plt.subplots(nrows=2, ncols=1, figsize=(6, 5), sharex=True,
                                 gridspec_kw={'height_ratios': [6, 1]})
        ax4[0].plot(f, r, 'k.', markersize=12)

    rmin = num.min(r)
    rmax = num.max(r)
    fac_norm = rmin/df

    flat_areas = []

    for i in range(int(num.floor(n/2)), int(len(r)-num.floor(n/2))):

        if plot_flat_range is True:
            reg, p, xp = calc_reg_prepplot(n, i, f, r, deg,
                                           plot=plot_flat_range)
            ax4[0].plot(xp, p(xp), '--')

        else:
            reg = calc_reg_prepplot(n, i, f, r, deg, plot=plot_flat_range)
        if abs(reg[0])/fac_norm < 1:
            flat_areas.append(i)

    _flat_f_ranges = get_flat_freq_ranges(r, f, n, flat_areas)

    if plot_flat_range:
        return _flat_f_ranges, rmin, rmax, fig4, ax4
    else:
        return _flat_f_ranges, rmin, rmax


def plot_f_range(flat_f_ranges, yval, fig, ax):
    for fl in flat_f_ranges:
        l_fr = num.linspace(fl[0], fl[1], 100)
        l_y = [1 for i in range(len(l_fr))]
        ax[1].plot(l_fr, l_y, 'r.', markersize=2)
        ax[1].yaxis.set_ticks([])
    return fig, ax


class dict_stats(Object):
    '''
    Dict for all stations + their flat freq ranges
    '''
    FlatFreqRanges = Dict.T(String.T(), List.T(Tuple.T(2, Float.T())))


def dump_flat_ranges(flat_f_ranges_stlist, nslc_list):
    '''
    writing flat areas of all stations to file

    '''
    fr_dict = dict(zip(map(lambda x: x,
                           nslc_list),
                           flat_f_ranges_stlist))
    flfr = dict_stats(FlatFreqRanges=fr_dict)
    flfr.regularize()
    flfr.validate()

    flfr.dump(filename='testing_psd_flat_ratio_alparray.yaml')


def const_psd_rat(mean_rat, cha, st, f_syn_keep,
                  dir_f=False, plot_flat_range=False):
    '''
    Find frequency range of roughly constant relative gain factor
    (optimal is no needed gain factor) in mean psd ratio of syn and
    obs data.

    :param mean_rat: Mean of PSD ratios for current nslc over all events
    '''

    net = st.network
    stat = st.station
    # num.save('mn_rats_%s_%s_%s.np' % (net, stat, cha), mean_rat)
    # num.save('freqs', f_syn_keep)

    # data
    freqs = f_syn_keep
    rat_a = mean_rat
    assert len(freqs) == len(rat_a)

    # settings for poloynomial fit
    n = 35
    deg = 1

    if plot_flat_range:
        _flat_f_ranges, rmin, rmax, fig4, ax4 = get_flat_freqs(n, freqs,
                                                               rat_a, deg=deg,
                                                plot_flat_range=plot_flat_range)
    else:
        _flat_f_ranges, rmin, rmax = get_flat_freqs(n, freqs, rat_a, deg=1)

    flat_f_ranges = [i for i in _flat_f_ranges if i[0] != i[1]]
    # print(flat_f_ranges)
    if plot_flat_range is True:
        yval = -100  # rmax + rmax/10
        fig4, ax4 = plot_f_range(flat_f_ranges, yval, fig4, ax4)
        fig4.savefig('%s/testing/%s_%s_%s_flatrange.png'
                     % (dir_f, net, stat, cha))
        # plt.show(fig4)
        plt.close(fig4)

    return flat_f_ranges


def plot_mean_ratio(mean_rat, f_syn_Z_nonan, nsl, dir_f, cha):
    '''
    Plotting of mean ratio of syn and obs PSD over all events
    (for current channel).

    '''
    fig3, ax3 = plt.subplots(nrows=1, ncols=1, figsize=(6, 5))

    if not mean_rat == []:
        ax3.plot(f_syn_Z_nonan[1:], mean_rat[1:], 'k.')
        fmax = max(f_syn_Z_nonan)
        fmin = 0.001
        ax3.set_xlim(fmin, fmax)
        ax3.set_xscale('log')
        # ax3.set_yscale('log')
        fig3.savefig('%s/testing/%s_%s_%s_meanratio_3.png'
                     % (dir_f, str(nsl[0]), str(nsl[1]), cha))
        plt.close(fig3)


def calc_mean_ratio(ratio_npar):   # hier dran arbeiten!
    '''
    Calculate mean of ratios for one station over all events.
    '''
    mean_rat = num.nanmean(ratio_npar, axis=0)
    return(mean_rat)


def plot_ratio_psd(fig2, ax2, i_x, i_y, nrows,
                   rats, freqs, ev_time_str):
    '''
    Plotting of ratio between obs and syn PSD for single events.
    per default not called.
    '''

    ax2[i_x, i_y].set_title(ev_time_str)
    ax2[i_x, i_y].plot(freqs[1:], rats[1:], 'b.')
    ymin, ymax = ax2[i_x, i_y].get_ylim()

    if i_x == nrows-1:
        ax2[i_x, i_y].set_xlabel('Frequency [Hz]')

    if i_y == 0:
        ax2[i_x, i_y].set_ylabel('PSD ratio (syn & obs)')

    fmax = max(max(freqs[1:]), max(freqs[1:]))
    fmin = 0.001
    ax2[i_x, i_y].set_xlim(fmin, fmax)
    ax2[i_x, i_y].set_ylim(ymin, ymax)
    ax2[i_x, i_y].set_xscale('log')


def plot_obs_syn_psd(ax, i_x, i_y, nrows, ev_time_str,
                     f_obs_Z, a_obs_Z, f_syn_Z, a_syn_Z,
                     nplots):
    '''
    Plotting of PSD for single events.

    '''

    ax[i_x, i_y].set_title(ev_time_str)
    ax[i_x, i_y].set_xscale('log')
    # ax[i_x, i_y].set_yscale('log')

    if i_x == nrows-1:
        ax[i_x, i_y].set_xlabel('Frequency [Hz]')

    if i_y == 0:
        ax[i_x, i_y].set_ylabel('PSD')

    if len(f_obs_Z) > 2 and len(a_obs_Z) > 2:
        try:
            ax[i_x, i_y].plot(f_obs_Z[1:], a_obs_Z[1:], 'k')
            nplots += 1

        except ValueError:
            print('ValueError obs')

    if len(f_syn_Z) > 2 and len(a_syn_Z) > 2:
        try:
            ax[i_x, i_y].plot(f_syn_Z[1:], a_syn_Z[1:], 'r')
            nplots += 1

        except ValueError:
            print('ValueError syn')

    if len(f_obs_Z) > 2 and len(a_obs_Z) > 2\
      and len(f_syn_Z) > 2 and len(a_syn_Z) > 2:

        # amin = min(min(a_obs_Z[1:]), min(a_syn_Z[1:]))
        # amax = max(max(a_obs_Z[1:]), max(a_syn_Z[1:]))
        # ax[i_x, i_y].set_ylim(amin, amax)
        fmax = max(max(f_syn_Z[1:]), max(f_obs_Z[1:]))
        fmin = 0.001
        ax[i_x, i_y].set_xlim(fmin, fmax)

    return nplots


def get_a_f(traces, cha):
    '''
    Calculation of actual PSD (freqs and ampls).

    :returns: 1-D numpy arrays f and a with frequencies and amplitudes.
    '''

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
        a[0] /= 2.  # erster wert durch zwei geteilt. warum?
        a[a.size//2] /= 2.  # mittlerer wert durch zwei geteilt. warum?
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
                   dir_f, arrT_array, st, i_st,
                   plot_psds=False, plot_ratio_extra=False,
                   plot_mean_rat=False, plot_flat_ranges=False):
    '''
    Next level function, called by ```prep_psd_fct``` to procede with current
    synthetic and observed datapiles.

    Selects for each channel and event traces from the piles and passes them
    on to function ```get_a_f``` that calculates the PSD and to plotting
    functions.

    For parameters see ```prep_psd_fct```.

    '''
    n_ev = len(catalog)
    nrows = math.ceil(n_ev / 5)
    ncols = 5

    freq_rat_list = []
    nslc_list = []

    for cha in ['Z', 'R', 'T']:
        if plot_psds is True:
            fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(13, 9))
            nplots = 0
        if plot_ratio_extra:
            fig2, ax2 = plt.subplots(nrows=nrows, ncols=ncols, figsize=(13, 9))

        cnt = 0

        for i_ev, ev in enumerate(catalog):
            ev_time_str = util.time_to_str(ev.time)[0:10]
            # print(ev_time_str)
            # abs. arrival times: arrT_array,
            # shape (len(subset_catalog), len(ns))

            start_twd = arrT_array[i_ev, i_st] - 60  # s
            end_twd = start_twd + 1800

            nsl = st.nsl()
            trs_obs = data_pile.all(
                tmin=start_twd,
                tmax=end_twd,
                tinc=400,                   # ???
                tpad=200,                    # ???
                trace_selector=lambda tr: tr.nslc_id[:2] == nsl[:2] and\
                                          tr.nslc_id[3] == cha,
                want_incomplete=True)

            trs_syn = syn_data_pile.all(
                tmin=start_twd,
                tmax=end_twd,
                tinc=400,                   # ???
                tpad=200,                    # ???
                trace_selector=lambda tr: tr.nslc_id[:2] == nsl[:2] and\
                                          tr.nslc_id[3] == cha,
                want_incomplete=True)

            if trs_syn and trs_obs:
                f_obs_Z, a_obs_Z = get_a_f(trs_obs, cha)
                f_syn_Z, a_syn_Z = get_a_f(trs_syn, cha)

                i_x = int(i_ev/ncols)
                i_y = int(i_ev % ncols)

                if plot_psds is True:
                    nplots = plot_obs_syn_psd(ax, i_x, i_y, nrows,
                                              ev_time_str,
                                              f_obs_Z, a_obs_Z,
                                              f_syn_Z, a_syn_Z,
                                              nplots)

                if num.array_equal(f_syn_Z, f_obs_Z) and len(a_syn_Z) > 0:
                    # diff_a = abs(a_syn_Z - a_obs_Z)

                    ratio_a = [a1/a2 if not a2 == 0 else 0
                               for a1, a2 in zip(a_syn_Z, a_obs_Z)]

                    # print([(f, a1, a2, a1/a2, d, '\n') for f,a1,a2,d
                    #        in zip(f_obs_Z, a_obs_Z, a_syn_Z, ratio_a)])
                    if plot_ratio_extra:
                        plot_ratio_psd(fig2, ax2, i_x, i_y, nrows,
                                       ratio_a, f_syn_Z,
                                       ev_time_str)
                    # ax[i_x, i_y].plot(f_syn_Z[1:], ratio_a[1:], 'b.')
                    if cnt == 0:
                        ratio_npar = num.empty([n_ev, len(f_syn_Z)])
                        ratio_npar[:] = num.nan

                    ratio_npar[i_ev, :] = ratio_a
                    cnt += 1
                    f_syn_keep = f_syn_Z
        if plot_psds:
            fig.tight_layout()
        if plot_ratio_extra:
            fig2.tight_layout()
        if plot_psds is True:
            if nplots != 0:
                fig.savefig('%s/testing/%s_%s_%s_3.png'
                            % (dir_f, str(nsl[0]), str(nsl[1]), cha))
            if plot_ratio_extra:
                fig2.savefig('%s/testing/%s_%s_%s_ratio_3.png'
                             % (dir_f, str(nsl[0]), str(nsl[1]), cha))

            plt.close(fig)
            if plot_ratio_extra:
                plt.close(fig2)

        if ratio_npar.shape != (1,):
            mean_rat = calc_mean_ratio(ratio_npar)
            # print(mean_rat)
            if plot_mean_rat:
                plot_mean_ratio(mean_rat, f_syn_keep, nsl, dir_f, cha)

            freq_rat_list.append(const_psd_rat(mean_rat, cha, st, f_syn_keep,
                                               plot_flat_range=plot_flat_ranges,
                                               dir_f=dir_f))
            nslc_list.append('%s.%s.%s.%s' % (nsl[0], nsl[1], nsl[2], cha))

    return freq_rat_list, nslc_list


def prep_psd_fct(i_st, st, subset_catalog, dir_f, arrT_array, datapath,
                 syndatapath, plot_psds=False, plot_ratio_extra=False,
                 plot_mean_rat=False, plot_flat_ranges=False):
    '''
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
    '''

    st_name = st.station
    st_data_pile = pile.make_pile(datapath, regex='_%s_' % st_name,
                                  show_progress=False)
    st_syn_data_pile = pile.make_pile(syndatapath, regex='_%s_' % st_name,
                                      show_progress=False)
    freq_rat_list = []
    nslc_list = []

    if st_data_pile.tmin is not None and st_data_pile.tmax is not None and\
      st_syn_data_pile.tmin is not None\
        and st_syn_data_pile.tmax is not None:

            freq_rat_list, nslc_list = calc_plot_psds(subset_catalog,
                                                      st_data_pile,
                                                      st_syn_data_pile,
                                                      dir_f, arrT_array,
                                                      st, i_st,
                                                      plot_psds=plot_psds,
                                                      plot_ratio_extra=plot_ratio_extra,
                                                      plot_mean_rat=plot_mean_ratio,
                                                      plot_flat_ranges=plot_flat_ranges)

    return freq_rat_list, nslc_list


# open questions --> sebastian: for some events offset between syn
# and obs psd,
# wrong mannitude?
# settings of tr.transfer etc
