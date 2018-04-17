import numpy as num
import math

import matplotlib.pyplot as plt
from pyrocko import util, trace


def get_a_f(traces):
    a_list = []

    for tr in traces:
        if tr.channel not in ('HHZ'):  # erstmal nur Z! vorher sonst rotieren!
            continue
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
        print('ValueError')
        f = []
        a = []
        return f, a


def calc_plot_psds(stat_list, catalog, data_pile,
                   syn_data_pile, dir_f, arrT_array):
    n_ev = len(catalog)
    nrows = math.ceil(n_ev / 5)
    ncols = 5

    for i_st, st in enumerate(stat_list):
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(13, 9))
        nplots = 0

        for i_ev, ev in enumerate(catalog):

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
                trace_selector=lambda tr: tr.nslc_id[:2] == nsl[:2],
                want_incomplete=True)
            trs_syn = syn_data_pile.all(
                tmin=start_twd,
                tmax=end_twd,
                tinc=400,                   ### ???
                tpad=200,                    ### ???
                trace_selector=lambda tr: tr.nslc_id[:2] == nsl[:2],
                want_incomplete=True)

            # trace.snuffle(trs_obs)
            # trace.snuffle(trs_syn)
            # trace.snuffle(trs_obs + trs_syn)

            f_obs_Z, a_obs_Z = get_a_f(trs_obs)
            f_syn_Z, a_syn_Z = get_a_f(trs_syn)

            i_x = int(i_ev/ncols)
            i_y = int(i_ev % ncols)

            ax[i_x, i_y].set_title(util.time_to_str(ev.time)[0:10])
            ax[i_x, i_y].set_xscale('log')
            ax[i_x, i_y].set_yscale('log')

            if i_x == nrows-1:
                print(i_x)
                ax[i_x, i_y].set_xlabel('Frequency [Hz]')

            if i_y == 0:
                ax[i_x, i_y].set_ylabel('PSD')

            if len(f_obs_Z) > 2 and len(a_obs_Z) > 2:
                try:
                    ax[i_x, i_y].plot(f_obs_Z[1:], a_obs_Z[1:], 'k')
                    nplots += 1

                except ValueError:
                    print('ValueError')

            if len(f_syn_Z) > 2 and len(a_syn_Z) > 2:
                try:
                    ax[i_x, i_y].plot(f_syn_Z[1:], a_syn_Z[1:], 'r')
                    nplots += 1

                except ValueError:
                    print('ValueError')

            if len(f_obs_Z) > 2 and len(a_obs_Z) > 2\
            and len(f_syn_Z) > 2 and len(a_syn_Z) > 2:

                amin = min(min(a_obs_Z[1:]), min(a_syn_Z[1:]))
                amax = max(max(a_obs_Z[1:]), max(a_syn_Z[1:]))
                ax[i_x, i_y].set_ylim(amin, amax)
                fmax = max(max(f_syn_Z[1:]), max(f_obs_Z[1:]))
                ax[i_x, i_y].set_xlim(0.01, fmax)

        plt.tight_layout()

        #if nplots != 0:
        #    plt.savefig('%s/%s_%s.png' % (dir_f, str(nsl[0]), str(nsl[1])))

        plt.close()
