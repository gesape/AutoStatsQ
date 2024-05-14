FROM ubuntu:20.04

MAINTAINER Gesa Petersen

# Update the package index
RUN apt-get -y update && apt-get upgrade -y -o Dpkg::Options::="--force-confold"

# Install required dependencies
RUN apt install -y build-essential
RUN DEBIAN_FRONTEND=noninteractive apt install -y cmake
RUN apt install -y libcurl4-gnutls-dev
RUN apt install -y libnetcdf-dev

# Install optional dependencies
RUN apt install -y gdal-bin libgdal-dev libfftw3-dev libpcre3-dev liblapack-dev libblas-dev libglib2.0-dev ghostscript

# to enable movie-making
RUN apt install -y graphicsmagick ffmpeg poppler-utils graphicsmagick-imagemagick-compat

# to enable document viewing via gmt docs
RUN apt install -y xdg-utils

# to build the documentation
RUN apt install -y python3-sphinx nano vim

# Install extras for Pyrocko
RUN apt install -y wget git

# Install GMT
# RUN cd /root && wget "https://github.com/GenericMappingTools/gmt/releases/download/5.4.5/gmt-5.4.5-src.tar.gz" && gunzip gmt-5.4.5-src.tar.gz && tar -xf gmt-5.4.5-src.tar
# RUN cd /root/gmt-5.4.5 && mkdir build
# RUN cd /root/gmt-5.4.5/build && cd /root/gmt-5.4.5/build && cmake .. && cmake --build . && cmake --build . --target install
RUN apt install -y gmt

# Install Pyrocko
# RUN cd /root && git clone https://git.pyrocko.org/pyrocko/pyrocko.git pyrocko
# RUN apt-get install -y make git python3-dev python3-setuptools python3-pip python3-wheel python3-numpy python3-numpy-dev python3-scipy \
#     python3-matplotlib python3-pyqt5 python3-pyqt5.qtopengl python3-pyqt5.qtsvg python3-pyqt5.qtwebkit python3-yaml \
#     python3-progressbar python3-jinja2 python3-requests python3-coverage python3-nose
# RUN cd /root/pyrocko && python3 -m pip install --no-deps --no-build-isolation --force-reinstall .
RUN apt -y install python3-pip
RUN pip3 install --upgrade pip
RUN pip3 install pyrocko[gui]

# Install Grond
RUN pip3 install grond

# Install GSHH and link to proper path
RUN apt install -y gmt-gshhg && ln -s /usr/share/gmt-gshhg /usr/local/share/coast

# Add sysop user
RUN /usr/sbin/useradd -m -s /bin/bash -p '*' sysop

# Install AutoStatsQ
RUN cd /home/sysop && su - sysop -c "git clone https://github.com/gesape/AutoStatsQ" && cd AutoStatsQ && python3 setup.py install

# Download Green functions (optional to compare amplitudes to synthetic amplitudes)
# RUN cd /home/sysop && su sysop -c "fomosto download kinherd global_2s"
