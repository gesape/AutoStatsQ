from pyrocko import pile, trace, util
import numpy as num
import math
import matplotlib.pyplot as plt

def test(input):
    print(input)


def ccs_allstats_one_event(i_ev, ev, stat_list, datapath, syndatapath,
                           out_dir, bp, arrT_array, cc_thresh,
                           debug_mode=False):
    '''
    for one event: call cc_single_stat_single_event for each station,
    collect optimal time shifts
    return list with timeshift, fixed order of stations!
    '''

    ev_t_str = util.time_to_str(ev.time).replace(' ', '_')
    # load data and syndata
    p_obs = pile.make_pile(datapath+ev_t_str, show_progress=False)
    p_syn = pile.make_pile(syndatapath+ev_t_str, show_progress=False)

    tshift_list = []

    if p_obs and p_syn:
        for i_st, st in enumerate(stat_list):
            # print(st.network, st.station)
            tmin = arrT_array[i_ev, i_st]-30
            tr_obs = p_obs.all(
                            trace_selector=lambda tr: tr.network == st.network
                            and tr.station == st.station
                            and tr.location not in ['10', '60']
                            and tr.channel == 'Z',
                            tmin=tmin,
                            tmax=tmin+300,                            
                            want_incomplete=True)
            tr_syn = p_syn.all(
                            trace_selector=lambda tr: tr.network == st.network
                            and tr.station == st.station
                            and tr.channel == 'Z',
                            tmin=tmin,
                            tmax=tmin+300,
                            want_incomplete=True)
            if len(tr_obs) != 0 and len(tr_syn) != 0:
                tr_syn = tr_syn[0]
                tr_obs = tr_obs[0]
                tr_obs.bandpass(bp[0], bp[1], bp[2])
                tr_syn.bandpass(bp[0], bp[1], bp[2])

                c = trace.correlate(tr_syn, tr_obs, mode='same', normalization='normal')
                t, coef = c.max()

                if debug_mode is True:
                    print(t, coef)
                    trace.snuffle([tr_syn, tr_obs])
                    trace.snuffle([c])

                if coef > cc_thresh:
                    tshift_list.append(t)
                else:
                    tshift_list.append(num.nan)
            else:
                tshift_list.append(num.nan)

    return tshift_list

def correct_for_med_tshifts(tshift_array):
    '''
    get median time shift of each event
    subtract median of each single value
    return corrected array
    '''
    tshift_medians = num.nanmedian(tshift_array, axis=0)
    print(tshift_medians)
    print(tshift_medians.shape)

    tshift_corr = num.empty((tshift_array.shape))
    print(tshift_corr)
    print(tshift_corr.shape)

    for i_ev in range(tshift_medians.shape[0]):
        for i_st in range(tshift_array.shape[0]):
            v = tshift_array[i_st, i_ev]
            if not num.isnan(v):
                tshift_corr[i_st, i_ev] = v - tshift_medians[i_ev]
            else:
                tshift_corr[i_st, i_ev] = v

    return tshift_corr



def plot_tshifts(tshifts_cor, means, stdevs, outfile, stations):
    # scatter plot: for each station all offsets
    n_plotrows = 1#math.ceil(len(stations)/20)
    fig, ax = plt.subplots(nrows=1, figsize=(20, n_plotrows*10))


    #for i_row in range(n_plotrows):
    ax.set_xlim(0, len(stations))
    ax.set_ylim(num.nanmin(tshifts_cor), num.nanmax(tshifts_cor))
    for i_st, st in enumerate(stations):
        yval = tshifts_cor[i_st, :]
        xval = [i_st for val in yval]
        # print(yval, xval)
        ax.scatter(xval, yval, marker='o')

            #    if not num.isnan(tshifts_cor[i_st, i_ev]):
            #        print(i_st, i_ev, tshifts_cor[i_st, i_ev])
            #        ax[i_row].plot(tshifts_cor[i_st, i_ev], [i_st])
        ax.plot(i_st, means[i_st], 'r.')


    ax.set_ylabel('Time shift [s]')
    stats = ['%s.%s' % (st.network, st.station) for st in stations]
    ax.set_xticks(num.arange(len(stats)))
    ax.set_xticklabels(stats, rotation=60)
    plt.tight_layout()
    fig.savefig(outfile)
    plt.close()
