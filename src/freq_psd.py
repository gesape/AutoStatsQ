import numpy as num
import math

import matplotlib.pyplot as plt
from pyrocko import util, trace, pile


def get_a_f(traces, cha):
    a_list = []

    for tr in traces:
        #print(tr.channel)
        #if tr.channel not in ('Z', 'R','T'):  # erstmal nur Z! vorher sonst rotieren!
        #    continue
        win = num.ones(tr.data_len())
        #tr.snuffle()

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
        #print('ValueError %s' %cha)
        f = []
        a = []
        return f, a

def plot_mean_ratio(ratio_npar, f_syn_Z, nsl, dir_f, cha):
    fig3, ax3 = plt.subplots(nrows=1, ncols=1, figsize=(6, 5))
    sumrat=0
    #print(nsl)
    #for i in range(19):
    #    print(ratio_npar[i,1])
    #    sumrat += ratio_npar[i,1]
    #print('hand mean', sumrat/19)
    mean_rat = num.mean(ratio_npar, axis=0)
    print('np mean', mean_rat[1])
    #print(mean_rat.shape)
    f_syn_Z_nonan = []
    mean_rat_nonan = []
    for f, m in zip(f_syn_Z[1:], mean_rat[1:]):
        if m > 0:
            f_syn_Z_nonan.append(f)
            mean_rat_nonan.append(m)
    if not mean_rat == []:
        ax3.plot(f_syn_Z_nonan, mean_rat_nonan, 'k.')
        fmax = max(f_syn_Z_nonan)
        fmin = 0.01
        ax3.set_xlim(fmin, fmax)
        ax3.set_xscale('log')    
        fig3.savefig('%s/testing/%s_%s_%s_meanratio.png' % (dir_f, str(nsl[0]), str(nsl[1]), cha))
        plt.close(fig3)

def plot_ratio_psd(fig2, ax2, i_x, i_y, nrows,
                   rats, freqs, ev_time_str):

    ax2[i_x, i_y].set_title(ev_time_str)
    ax2[i_x, i_y].plot(freqs[1:], rats[1:], 'b.')
    ymin, ymax = ax2[i_x, i_y].get_ylim()

    if i_x == nrows-1:
        ax2[i_x, i_y].set_xlabel('Frequency [Hz]')

    if i_y == 0:
        ax2[i_x, i_y].set_ylabel('PSD ratio (syn & obs')

    fmax = max(max(freqs[1:]), max(freqs[1:]))
    fmin = 0.01
    ax2[i_x, i_y].set_xlim(fmin, fmax)
    ax2[i_x, i_y].set_ylim(ymin, ymax)
    ax2[i_x, i_y].set_xscale('log')



def plot_obs_syn_psd(ax, i_x, i_y, nrows, ev_time_str,
                     f_obs_Z, a_obs_Z, f_syn_Z, a_syn_Z,
                     nplots):

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
            #print(len(f_syn_Z), len(a_syn_Z))
            ax[i_x, i_y].plot(f_syn_Z[1:], a_syn_Z[1:], 'r')
            nplots += 1

        except ValueError:
            print('ValueError syn')

    if len(f_obs_Z) > 2 and len(a_obs_Z) > 2\
      and len(f_syn_Z) > 2 and len(a_syn_Z) > 2:

        amin = min(min(a_obs_Z[1:]), min(a_syn_Z[1:]))
        amax = max(max(a_obs_Z[1:]), max(a_syn_Z[1:]))
        ax[i_x, i_y].set_ylim(amin, amax)
        fmax = max(max(f_syn_Z[1:]), max(f_obs_Z[1:]))
        fmin = 0.01
        ax[i_x, i_y].set_xlim(fmin, fmax)

    return nplots


def calc_plot_psds(catalog, data_pile, syn_data_pile,
                   dir_f, arrT_array, st, i_st):
    n_ev = len(catalog)
    nrows = math.ceil(n_ev / 5)
    ncols = 5

    for cha in ['Z','R','T']:
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(13, 9))
        fig2, ax2 = plt.subplots(nrows=nrows, ncols=ncols, figsize=(13, 9))

        nplots = 0
        # print(st.nsl)
        cnt = 0
        ratio_npar = 0
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
                trace_selector=lambda tr: tr.nslc_id[:2] == nsl[:2] and tr.nslc_id[3] == cha,
                want_incomplete=True)

            trs_syn = syn_data_pile.all(
                tmin=start_twd,
                tmax=end_twd,
                tinc=400,                   ### ???
                tpad=200,                    ### ???
                trace_selector=lambda tr: tr.nslc_id[:2] == nsl[:2] and tr.nslc_id[3] == cha,
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

                if len(f_syn_Z) == len(f_obs_Z):
                    diff_a = abs(a_syn_Z - a_obs_Z)
                    ratio_a = abs(a_syn_Z / a_obs_Z)
                    #print([(f, a1, a2, d) for f,a1,a2,d
                    #in zip(f_obs_Z, a_obs_Z, a_syn_Z, diff_a,)])
                    plot_ratio_psd(fig2, ax2, i_x, i_y, nrows,
                                   ratio_a, f_syn_Z,
                                   ev_time_str)

                    if cnt == 0:
                        ratio_npar = num.empty([n_ev, len(f_syn_Z)])

                    ratio_npar[i_ev]=ratio_a
                    cnt += 1
                    f_syn_keep = f_syn_Z

                
        fig.tight_layout()
        fig2.tight_layout()

        if nplots != 0:
            fig.savefig('%s/testing/%s_%s_%s_t.png' % (dir_f, str(nsl[0]), str(nsl[1]), cha))
            fig2.savefig('%s/testing/%s_%s_%s_ratio.png' % (dir_f, str(nsl[0]), str(nsl[1]), cha))

        plt.close(fig)
        plt.close(fig2)

        if ratio_npar.shape != (1,):

            plot_mean_ratio(ratio_npar, f_syn_keep, nsl, dir_f, cha)

        

def prep_psd_fct(i_st, st, subset_catalog, dir_f, arrT_array, datapath, syndatapath):
    st_name = st.station
    st_data_pile = pile.make_pile(datapath, regex='_%s_' %st_name,
                                  show_progress=False)
    st_syn_data_pile = pile.make_pile(syndatapath, regex='_%s_' %st_name,
                                      show_progress=False)

    if st_data_pile.tmin is not None and st_data_pile.tmax is not None and\
          st_syn_data_pile.tmin is not None and st_syn_data_pile.tmax is not None:

            calc_plot_psds(subset_catalog,
                           st_data_pile,
                           st_syn_data_pile,
                           dir_f, arrT_array,
                           st, i_st)



# calc some trend or averaged or smoothed curve for syn and obs psd, 
# calc. difference
# plot difference
# determine f range, for which gain factor is constant