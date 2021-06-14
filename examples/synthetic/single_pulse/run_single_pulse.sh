# File: run_single_pulse.sh

# This synthetic dataset contains the following

#  - 16 Seismic Analysis Code (SAC) files that contain the waveforms
#  -  1 stationtable file of infrasound array DIA (see Evers, 2008)

# ------------

# The following signals are contained in the data:

# 1. A simple Ricket wavelet from 90 deg (E), propagating with 340 m/s
#    This signal has a peak amplitude of 10 pascal and is mostly coherent between 1.0-5.0 Hz.

# 2. Incoherent noise with amplitude 1.0 pascal

# ------------
# Using the following script, process the full recording with 20 second time windows with 90% overlap:

script=../../../python/array_process.py
array_name='NL.DIA'
inventory='../../NL.DIA.xml'
deconvolution='gain' # or 'response'
skip_elements=''

starttime='2007-01-01T00:00:00'
endtime='2007-01-01T00:30:00'

tbinsize='20'
overlap='90'
bandpass='0.1/5.0'

data='data'
srate='40'

azimuth_range='0./359./1.0'
appvel_range='300.0/400.0/5.0'
grid_type='tele'
grid_return='max'

#Time-domain Fisher
${script} \
    -n ${array_name} -m timefisher \
    -st ${starttime} -et ${endtime} \
    -bs ${tbinsize} -ov ${overlap} \
    -fr ${bandpass} -inv ${inventory} \
    -th ${azimuth_range} -ct ${appvel_range} \
    -gt ${grid_type} -gr ${grid_return} \
    -dr ${data} -fs ${srate} -el ${skip_elements} \
    -decon ${deconvolution}

# You can further enhance detectability of signal 2 by applying a frequency filter and/or by
# 'beaming' in the direction of the signal of interest (this corresponds to a wavenumber filter).
# For example, it is possible to limit the back azimuth angles to 45-135 degrees:

# azimuth_range='45./135./1.0'

# #Time-domain Fisher
# ${script} \
#     -n ${array_name} -m timefisher \
#     -st ${starttime} -et ${endtime} \
#     -bs ${tbinsize} -ov ${overlap} \
#     -fr ${bandpass} -inv ${inventory} \
#     -th ${azimuth_range} -ct ${appvel_range} \
#     -gt ${grid_type} -gr ${grid_return} \
#     -dr ${data} -fs ${srate} -el ${skip_elements} \
#     -decon ${deconvolution}

# Next step is Fisher FK analysis, to analyse the coherence as a function of time and frequency.
# The following script analyses in windows of 30 seconds with 50% overlap, over a band from
# 0.05-10 Hz, with spectral steps of 0.005 Hz. The results are averaged over 10 frequency bands.

fbinsize='30'
overlap='50'
spectral_band='0.05/10.0/0.005/10'

# FK+ Fisher (script uses beamforming results from timefisher)
${script} \
    -n ${array_name} -m tfreqfisher \
    -st ${starttime} -et ${endtime} \
    -bs ${fbinsize} -ov ${overlap} \
    -fr ${spectral_band} -inv ${inventory} \
    -dr ${data} -fs ${srate} -el ${skip_elements} \
    -decon ${deconvolution} -tdm timefisher

# # FK Fisher (script beamforms over grid specified above)
# ./array_process.py \
#     -n ${array_name} -m freqfisher \
#     -st ${starttime} -et ${endtime} \
#     -bs ${fbinsize} -ov ${overlap} \
#     -fr ${spectral_band} -inv ${inventory} \
#     -th ${azimuth_range} -ct ${appvel_range} \
#     -gt ${grid_type} \
#     -dr ${data} -fs ${srate} -el ${skip_elements} \
#     -decon ${deconvolution}
