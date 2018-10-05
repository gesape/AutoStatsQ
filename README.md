# AutoStatsQ
Toolbox for automated station quality control for MT inversion (-work in progress-)
Please contact me for further description and help: gesap@gfz-potsdam.de

- Catalog search for teleseismic events with uniform azimuthal coverage around array
- Download of data & metadata for these events + computation of synthetic data
- Relative gain factors in time domain (relative to reference station)
- Rayleigh wave polarization analysis for detection of sensor misorientations
- Comparison of obs. and synth. PSDs; determining frequency ranges suitable for MT inversion


Requirements
------------

- python3
- seismology toolbox pyrocko: https://pyrocko.org/ (Heimann et al. 2017)

Download and Installation
-------------------------

cd into the folder where you want to do the installation\n
git clone https://github.com/gesape/AutoStatsQ \n
cd AutoStatsQ\n
(sudo) python3 setup.py install\n

Basic commands
--------------

show basic commands:


	autostatsq -h


generate an example config file:

	autostatsq --generate_config GENERATE_CONFIG

run 

	autostatsq --config name_of_config_file --run RUN

