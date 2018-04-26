from pyrocko import model, util, orthodrome


def subset_events_dist(catalog, mag_min, mag_max,
                       tmin, tmax, st_lat, st_lon,
                       dist_min=None, dist_max=None):
    '''
    Extract a subset of events from event catalog

    :param catalog: Event catalog in pyrocko format
    :param mag_min: Min. magnitude
    :param mag_max: Max. magnitude
    :param tmin: string representing UTC time
    :param tmax: string representing UTC time
    :param format tmin: time string format ('%Y-%m-%d %H:%M:%S.OPTFRAC')
    :param format tmax: time string format ('%Y-%m-%d %H:%M:%S.OPTFRAC')
    :param dist_min: min. distance (km)
    :param dist_max: max. distance (km)

    :returns: list of events
    '''

    use_events = []
    events = model.load_events(catalog)
    for ev in events:
        if ev.magnitude < mag_max and\
          ev.magnitude > mag_min and\
          ev.time < util.str_to_time(tmax) and\
          ev.time > util.str_to_time(tmin):
            if dist_min or dist_max:
                dist = orthodrome.distance_accurate50m_numpy(
                       ev.lat, ev.lon, st_lat, st_lon)/1000.

                if dist_min and dist_max and\
                  dist > dist_min and dist < dist_max:
                    use_events.append(ev)

                if dist_min and not dist_max and dist > dist_min:
                    use_events.append(ev)

                if dist_max and not dist_min and dist < dist_max:
                    use_events.append(ev)

    return(use_events)
