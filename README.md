# AutoStatsQ (-work in progress-)
Toolbox for automated station quality control for MT inversion 
Please contact me for further description and help: gesap@gfz-potsdam.de

- Catalog search for teleseismic events with uniform azimuthal coverage around array
- Download of data & metadata for these events + computation of synthetic data
- Relative gain factors in time domain (relative to reference station)
- Rayleigh wave polarization analysis for detection of sensor misorientations
- Comparison of obs. and synth. PSDs; determining frequency ranges suitable for MT inversion


Requirements
------------

- python3
- Seismology toolbox pyrocko: https://pyrocko.org/ (Heimann et al. 2017)
- To compute synthetic data a pre-calculated GF database can be downloaded using `fomosto download` from the pyrocko environment. http://kinherd.org:8080/gfws/static/stores/


Download and Installation
-------------------------

- cd into the folder where you want to do the installation
- git clone https://github.com/gesape/AutoStatsQ
- cd AutoStatsQ
- (sudo) python3 setup.py install

Basic commands
--------------

show basic commands/ help:

	autostatsq -h


generate an example config file:

	autostatsq --generate_config GENERATE_CONFIG

run 

	autostatsq --config name_of_config_file --run RUN


Config file settings
--------------------

The default config file should look like this (without the comments):

```yaml
--- !autostatsq.config.AutoStatsQConfig
Settings:
- !autostatsq.config.GeneralSettings
  data_dir: /some/data/directory/
  list_station_lists: [/path/to/station-file/file.csv, /path/to/station-file/file.xml]

- !autostatsq.config.CatalogConfig
  search_events: true 
  # search gCMT catalog for events?

  use_local_catalog: false
  # Or use a local (already downloaded) catalog? 
  # Needed re-runs using same catalog.

  subset_of_local_catalog: true
  # Find a subset of the full catalog?

  use_local_subsets: false
  # Use local (already saved) subset instead?
  
  subset_fns: {}
  # if so, give here paths to subset-catalog-files: e.g. 
  # {'deep': 'catalog_deep_subset.txt',
  # 'shallow': 'catalog_shallow_subset.txt'}

  catalog_fn: catalog.txt
  # filename of catalog that is downloaded or already exists

  min_mag: 6.5
  max_mag: 8.5
  tmin_str: '2000-01-01 00:00:00'
  tmax_str: '2018-10-01 00:00:00'
  min_dist_km: 1000.0
  max_dist_km: 9999999.9  
  depth_options:
    deep: [25000, 1000000] # [m]
    shallow: [100, 40000] # [m]

  wedges_width: 15
  # backazimuthal step for subset generation
  # adjust to get more events, especially if time range is small

  mid_point: [46.98, 10.74]
  # give a rough estimate of midpoint of array/ network

  ### catalog plotting options ###
  dist: 165.0
  # max. distance (degrees) for catalog plot  
  plot_catalog_all: false
  # plots entire catalog on a map

  plot_hist_wedges: false
  plot_wedges_vs_dist: false
  plot_wedges_vs_magn: false
  plot_dist_vs_magn: false
  # catalog statistics plot
  
  plot_catalog_subset: false
  # plots the subset(s) on a map

- !autostatsq.config.ArrTConfig
  calc_first_arr_t: true
  # Should first arrivals be computed?

  phase_select: P|p|P(cmb)P(icb)P(icb)p(cmb)p|P(cmb)Pv(icb)p(cmb)p|P(cmb)P<(icb)(cmb)p
  # which phases?

  calc_est_R: true
  # compute arrival time of Rayleigh waves for each station-event pair? 
  # (needed for orientation test)


- !autostatsq.config.CakeTTTGenerator
  # uses travel time tables instead of cake (faster, but settings more difficult)
  # if you are not familiar with fomosto' travel time tables use ArrTConfig section 
  # instead!
  calc_ttt: false
  dir_ttt: /directory/to/save/traveltimes/
  earthmodel_id: prem-no-ocean.f
  tabulated_phases:
  - !pf.TPDef
    id: p
    definition: P,p
  dist_min: 1000.0
  dist_max: 9999.0
  dist_acc: 2.0
  s_depth_min: 10.0
  s_depth_max: 50.0
  s_acc: 2.0
  r_depth_min: 0.0
  r_depth_max: 0.0
  r_acc: 2.0
  t_acc: 3.0

- !autostatsq.config.MetaDataDownloadConfig
  # download of metadata and data

  download_data: false
  download_metadata: false 
  use_downmeta: false
  # Set to true if downloaded metadata should be used.
  # If metadata with responses is stored locally, you can use  
  # local_metadata: [list_of_metadata_files] instead of or in 
  # addition to downloading
  components_download: HH*
  # '*' would download all and analyse the most broadband channel for each
  # station
  token:
    geofon: /home/gesap/Documents/AlpArray/download-waveforms/swathD/token.asc
  sites: [geofon, orfeus, iris]

  dt_start: 0.1
  # start time before origin time [h]
  dt_end: 1.5
  # end time after origin time [h]

- !autostatsq.config.RestDownRotConfig
  # restitution, downsampling and rotation of data
  # required for all tests
  rest_data: false
  freqlim: [0.005, 0.01, 0.2, 0.25] # [Hz]
  rotate_data: false
  deltat_down: 2

- !autostatsq.config.SynthDataConfig
  # computation of synthetic data
  # needed for PSD-test only, can otherwise be left out
  make_syn_data: false
  engine_path: /path/to/GF_stores
  store_id: global_2s

- !autostatsq.config.GainfactorsConfig
  # settings for first test
  calc_gainfactors: false
  gain_factor_method:
  - reference_nsl
  - [GE, MATE]
  ### describe different methods
  fband:
    corner_hp: 0.01 # [Hz]
    corner_lp: 0.2 # [Hz]
    order: 4
  taper_xfrac: 0.25 # [s]

  wdw_st_arr: 5
  wdw_sp_arr: 60
  # time window around P phase onset, start [s] before and end [s] after theo. arrival time

  phase_select: first(P|p|P(cmb)P(icb)P(icb)p(cmb)p|P(cmb)Pv(icb)p(cmb)p|P(cmb)P<(icb)(cmb)p)
  components: [Z, R, T]

  # plotting options
  plot_median_gain_on_map: false
  plot_allgains: false

- !autostatsq.config.PSDConfig
  # settings for PSD test
  calc_psd: false
  tinc: 600  # [s]
  tpad: 200  # [s]
  dt_start: 60  # [s]
  dt_end: 1800  # [s]
  n_poly: 25
  norm_factor: 50
  f_ign: 0.02  # [Hz]
  only_first: true
  plot_psds: false
  plot_ratio_extra: false
  plot_m_rat: false
  plot_flat_ranges: false
  plt_neigh_ranges: false

- !autostatsq.config.OrientConfig
  # settings for orientation test
  orient_rayl: false
  bandpass: [3.0, 0.01, 0.05]  # [Hz]
  start_before_ev: 30.0  # start befor theo. Rayleigh wave arrival, [s]
  stop_after_ev: 480.0  # end after theo. Rayleigh wave arrival, [s]
  ccmin: 0.80
  # min. cross-correlation value. results below this value will not be
  # considered
  plot_heatmap: false
  plot_distr: false
  plot_orient_map_fromfile: false

- !autostatsq.config.maps
  # settings for all output maps
  map_size: [30.0, 30.0]
  pl_opt: [46, 11.75, 800000]
  # mid point of map (lat, lon) and radius [m]
  pl_topo: false
  # plotting topography can be very slow,
  # topographic data will be downloaded first
```


Station list:
-------------

Lists of stations as input can be in station-xml format or as comma-spread-file with columns: network code, station code, latitude (float), longitude(float), station elevation [km], station depth [km].

Step-by-step instructions:
--------------------------

- make a working directory
- prepare station list for input
- generate a template config file (`autostatsq --generate_config GENERATE_CONFIG`) and adjust the settings
- it might be helpful to not run the toolbox in one run, but go though it step-by-step, setting the current step to `true` and the others to `false`... (`autostatsq --config my_config_file --run RUN`)

Output:
-------

Gain test:
- yaml files with median, mean and standard deviation of Ai,j/Ai,ref; i: event, j: station
- csv files with Ai,j/Ai,ref for all events and all stations
- optional: maps showing median log. Ai,j/Ai,ref for all stations and components
- optional: ...

Orientation test:
- yaml file with median, mean and standard deviation of obtained correction angle
- yaml file with polarity errors
- optional: map showing the median correction angle

PSD test:
- flat freq ranges in yaml file for each station and component
- optional: synth. and real PSDs
- optional: plots showing fit through PSD ratios 
