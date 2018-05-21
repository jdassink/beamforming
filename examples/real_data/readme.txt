# This real data dataset contains the following

#  - 13 Seismic Analysis Code (SAC) files that contain the waveforms
#  -  1 stationtable file of infrasound array DIA (see Evers, 2008)

# Note that 3 elements we not functioning and are not included in the stationtable either
# An analysis of this case can be found on page 44 of Evers, 2008:
# https://repository.tudelft.nl/islandora/object/uuid%3A4de38d6f-8f68-4706-bf34-4003d3dff0ce 

# To process the full recording with 12 second time windows with 90% overlap:

timefisher NL.DIA.dat 12 90 0 360 1 300 450 5 tele max sac *SAC
