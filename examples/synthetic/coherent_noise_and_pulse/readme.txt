# This synthetic dataset contains the following

#  - 16 Seismic Analysis Code (SAC) files that contain the waveforms
#  -  1 stationtable file of infrasound array DIA (see Evers, 2008)

------------

The following signals are contained in the data:

# 1. Coherent noise from 300 deg (NNW), propagating with 350 m/s
#    This signal has a signal noise ratio (SNR) of 1 and is coherent in the 0.2-0.3 Hz band

# 2. A simple Ricket wavelet from 90 deg (E), propagating with 340 m/s
#    This signal has a peak amplitude of 10 pascal and is mostly coherent between 1.0-5.0 Hz.

# 3. Incoherent noise, added so that signal 1 has a SNR of 1.

------------

# To process the full recording with 20 second time windows with 90% overlap:

timefisher NL.DIA.dat 20 90 0 360 1 300 400 5 tele max sac *SAC

# You can further enhance detectability of signal 2 by applying a frequency filter and/or by
# 'beaming' in the direction of the signal of interest (this corresponds to a wavenumber filter).
# For example, it is possible to limit the back azimuth angles to 45-135 degrees:

timefisher NL.DIA.dat 20 90 45 135 1 300 400 5 tele max sac *SAC
