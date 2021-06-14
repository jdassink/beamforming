array_name='IM.I32KE'
inventory='IM.I32KE.xml'
deconvolution='gain' # or 'response'
skip_elements=''

starttime='2021-05-22T16:00:00'
endtime='2021-05-23T00:00:00'
# starttime='2021-05-22T00:00:00'
# endtime='2021-05-25T00:00:00'

tbinsize='30'
overlap='50'
bandpass='0.4/0.6'

data='../data'
srate='20'

#Time-delay of arrival
# ./array_process.py \
#     -n ${array_name} -m tdoa \
#     -st ${starttime} -et ${endtime} \
#     -bs ${tbinsize} -ov ${overlap} \
#     -fr ${bandpass} -inv ${inventory} \
#     -dr ${data} -fs ${srate} -el ${skip_elements} \
#     -decon ${deconvolution}

azimuth_range='255./275./0.5'
appvel_range='320.0/360.0/2.0'
grid_type='tele'
grid_return='max'

#Time-domain Fisher
./array_process.py \
    -n ${array_name} -m timefisher \
    -st ${starttime} -et ${endtime} \
    -bs ${tbinsize} -ov ${overlap} \
    -fr ${bandpass} -inv ${inventory} \
    -th ${azimuth_range} -ct ${appvel_range} \
    -gt ${grid_type} -gr ${grid_return} \
    -dr ${data} -fs ${srate} -el ${skip_elements} \
    -decon ${deconvolution}

fbinsize='30'
overlap='50'
spectral_band='0.05/6.0/0.001/10'

# # FK Fisher
# ./array_process.py \
#     -n ${array_name} -m freqfisher \
#     -st ${starttime} -et ${endtime} \
#     -bs ${fbinsize} -ov ${overlap} \
#     -fr ${spectral_band} -inv ${inventory} \
#     -th ${azimuth_range} -ct ${appvel_range} \
#     -gt ${grid_type} \
#     -dr ${data} -fs ${srate} -el ${skip_elements} \
#     -decon ${deconvolution}

# FK+ Fisher (driven by timefisher)
./array_process.py \
    -n ${array_name} -m tfreqfisher \
    -st ${starttime} -et ${endtime} \
    -bs ${fbinsize} -ov ${overlap} \
    -fr ${spectral_band} -inv ${inventory} \
    -dr ${data} -fs ${srate} -el ${skip_elements} \
    -decon ${deconvolution} -tdm timefisher
