# File: run_real_data.sh

# This real data dataset contains the following

#  - 13 Seismic Analysis Code (SAC) files that contain the waveforms
#  -  1 stationtable file of infrasound array DIA (see Evers, 2008)

# Note that 3 elements we not functioning and are not included in the stationtable either
# An analysis of this case can be found on page 44 of Evers, 2008:
# https://repository.tudelft.nl/islandora/object/uuid%3A4de38d6f-8f68-4706-bf34-4003d3dff0ce 

# To process the full recording with 12 second time windows with 90% overlap:

script=../../python/array_process.py
array_name='NL.DIA'
inventory='../NL.DIA.xml'
deconvolution='' # or 'response'
skip_elements=''

starttime='2007-247T20:00:00'
endtime='2007-247T20:10:00'

tbinsize='12'
overlap='95'
bandpass='0.1/5.0'

data='data'
srate='40'

azimuth_range='0./359./0.5'
appvel_range='300.0/450.0/2.5'
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

# -------------------------------------------------------------------------------------------------
# Another method is the use of time-delay of arrival processing.
# The method relies on time-differences between sensor pairs to
# estimate the back azimuth and trace (apparent) velocity in a
# least-squares sense. This method works better with larger time windows

tbinsize='30'
overlap='95'
bandpass='0.1/5.0'

${script} \
    -n ${array_name} -m tdoa \
    -st ${starttime} -et ${endtime} \
    -bs ${tbinsize} -ov ${overlap} \
    -fr ${bandpass} -inv ${inventory} \
    -dr ${data} -fs ${srate} -el ${skip_elements} \
    -decon ${deconvolution}

# -------------------------------------------------------------------------------------------------

# Next step is Fisher FK analysis, to analyse the coherence as a function of time and frequency.
# The following script analyses in windows of 30 seconds with 50% overlap, over a band from
# 0.05-10 Hz, with spectral steps of 0.005 Hz. The results are averaged over 10 frequency bands.

fbinsize='30'
overlap='95'
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

