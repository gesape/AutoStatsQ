import numpy as num
import math

import matplotlib.pyplot as plt
from pyrocko import util, trace, pile


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


def plot_mean_ratio(ratio_npar, f_syn_Z, nsl, dir_f, cha):
    '''
    Plotting of mean ratio of syn and obs PSD over all events
    (for current channel).

    '''

    fig3, ax3 = plt.subplots(nrows=1, ncols=1, figsize=(6, 5))
    '''
    sumrat = 0
    dev = 0
    for i in range(19):
        if not num.isnan(ratio_npar[i,1]):
            print(ratio_npar[i,1])
            sumrat += ratio_npar[i,1]
            dev +=1
    if dev != 0:
        print('hand mean', sumrat/dev)
    '''
    mean_rat = num.nanmean(ratio_npar, axis=0)
    # print('np mean', mean_rat[1])
    f_syn_Z_nonan = []
    mean_rat_nonan = []
    for f, m in zip(f_syn_Z[1:], mean_rat[1:]):
        if m > 0:
            f_syn_Z_nonan.append(f)
            mean_rat_nonan.append(m)
    if not mean_rat == []:
        ax3.plot(f_syn_Z_nonan, mean_rat_nonan, 'k.')
        fmax = max(f_syn_Z_nonan)
        fmin = 0.001
        ax3.set_xlim(fmin, fmax)
        ax3.set_xscale('log')
        fig3.savefig('%s/testing/%s_%s_%s_meanratio.png'
                     % (dir_f, str(nsl[0]), str(nsl[1]), cha))
        plt.close(fig3)


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
    ax[i_x, i_y].set_yscale('log')

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


def calc_plot_psds(catalog, data_pile, syn_data_pile,
                   dir_f, arrT_array, st, i_st, plot_ratio_extra=False):
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

    for cha in ['Z', 'R', 'T']:
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(13, 9))
        if plot_ratio_extra:
            fig2, ax2 = plt.subplots(nrows=nrows, ncols=ncols, figsize=(13, 9))

        nplots = 0
        # print(st.nsl)
        cnt = 0
        ratio_npar = num.asarray((0,))  # just fir testing
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
                tinc=400,                   ### ???
                tpad=200,                    ### ???
                trace_selector=lambda tr: tr.nslc_id[:2] == nsl[:2] and\
                                          tr.nslc_id[3] == cha,
                want_incomplete=True)

            trs_syn = syn_data_pile.all(
                tmin=start_twd,
                tmax=end_twd,
                tinc=400,                   ### ???
                tpad=200,                    ### ???
                trace_selector=lambda tr: tr.nslc_id[:2] == nsl[:2] and\
                                          tr.nslc_id[3] == cha,
                want_incomplete=True)

            if trs_syn and trs_obs:
                f_obs_Z, a_obs_Z = get_a_f(trs_obs, cha)
                f_syn_Z, a_syn_Z = get_a_f(trs_syn, cha)

                i_x = int(i_ev/ncols)
                i_y = int(i_ev % ncols)

                nplots = plot_obs_syn_psd(ax, i_x, i_y, nrows,
                                          ev_time_str,
                                          f_obs_Z, a_obs_Z,
                                          f_syn_Z, a_syn_Z,
                                          nplots)

                if num.array_equal(f_syn_Z, f_obs_Z) and len(a_syn_Z) > 0:
                    #diff_a = abs(a_syn_Z - a_obs_Z)
                    ratio_a = [a1/a2 for a1, a2 in zip(a_syn_Z, a_obs_Z)]
                    #print([(f, a1, a2, a1/a2, d, '\n') for f,a1,a2,d
                    #        in zip(f_obs_Z, a_obs_Z, a_syn_Z, ratio_a)])
                    if plot_ratio_extra:
                        plot_ratio_psd(fig2, ax2, i_x, i_y, nrows,
                                       ratio_a, f_syn_Z,
                                       ev_time_str)
                    ax[i_x, i_y].plot(f_syn_Z[1:], ratio_a[1:], 'b.')
                    if cnt == 0:
                        ratio_npar = num.empty([n_ev, len(f_syn_Z)])
                        ratio_npar[:] = num.nan

                    ratio_npar[i_ev, :] = ratio_a
                    cnt += 1
                    f_syn_keep = f_syn_Z

        fig.tight_layout()
        if plot_ratio_extra:
            fig2.tight_layout()

        if nplots != 0:
            fig.savefig('%s/testing/%s_%s_%s.png'
                        % (dir_f, str(nsl[0]), str(nsl[1]), cha))
            if plot_ratio_extra:
                fig2.savefig('%s/testing/%s_%s_%s_ratio.png'
                             % (dir_f, str(nsl[0]), str(nsl[1]), cha))

        plt.close(fig)
        if plot_ratio_extra:
            plt.close(fig2)

        if ratio_npar.shape != (1,):

            plot_mean_ratio(ratio_npar, f_syn_keep, nsl, dir_f, cha)


def prep_psd_fct(i_st, st, subset_catalog, dir_f, arrT_array, datapath,
                 syndatapath, plot_ratio_extra=False):
    '''
    Preparing the PSD calculations and plotting, e.g. making data-piles
    and calling next function if both, synthetic and recorded data is
    in piles.

    :param i_st: number of station - must fit arrival time table!
    :param st: pyrocko station object
    :param subset_catalog: Catalog to be used (list of pyrocko events)
    :param dir_f: directory to store PSD results at
    :param arrT_array: Arrival times
    :param datapath: Path to restituted (,rotated and downsampled) data
    :param syndatapath: Path to sythetic data
    :param plot_ratio_extra: Should ratios of PSDs be plotted to extra figure?
                             They are displayed within the PSD plot anyway,
                             default is False.
    '''

    st_name = st.station
    st_data_pile = pile.make_pile(datapath, regex='_%s_' % st_name,
                                  show_progress=False)
    st_syn_data_pile = pile.make_pile(syndatapath, regex='_%s_' % st_name,
                                      show_progress=False)

    if st_data_pile.tmin is not None and st_data_pile.tmax is not None and\
          st_syn_data_pile.tmin is not None\
          and st_syn_data_pile.tmax is not None:

            calc_plot_psds(subset_catalog,
                           st_data_pile,
                           st_syn_data_pile,
                           dir_f, arrT_array,
                           st, i_st,
                           plot_ratio_extra=plot_ratio_extra)


# calc some trend or averaged or smoothed curve for syn and obs psd,
# calc. difference
# plot difference
# determine f range, for which gain factor is constant
