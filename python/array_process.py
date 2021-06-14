#!/usr/bin/env python
"""
Python driver program for various array processing routines

.. module:: array_process

:author:
    Jelle Assink (jelle.assink@knmi.nl)

:copyright:
    2021, Jelle Assink

:license:
    This code is distributed under the terms of the
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.en.html)
"""
import os
import sys
import shutil
import argparse
import numpy as np

from obspy import UTCDateTime, read_inventory, Stream
# from obspy.signal.filter import envelope
from wave_data_interface import get_data, select_inventory
from wave_data_interface import inv_to_df, df_to_ascii, get_gain

bf_methods = ['ccts', 'freqfisher', 'tdoa', 'timefisher', 'tfreqfisher']

def main(argv):
    print ('')
    print ('='*40)
    print (' --- Array Processing Wrapper Script ---')
    print ('='*40)
    print ('')
    # Get command line parameters
    proc = process_cli(argv)
    inv = read_inventory(proc.inventory)

    # Get waveform data
    data_source = dict(local=proc.data_root)
    st = get_data(data_source, inventory=inv,
                  starttime=proc.starttime,
                  endtime=proc.endtime)

    # Update inventory based on available data
    invs = select_inventory(st, inv,
                            starttime=proc.starttime,
                            endtime=proc.endtime)
    proc.__setattr__('inventory', invs)

    # Start processing
    t0 = proc.starttime
    t1 = proc.endtime
    ndays = (t1.datetime-t0.datetime).days

    if ndays >= 1:
        for day in range(0,ndays):
            dt0 = t0  + day*24*3600
            dt1 = dt0 + 24*3600
            proc.__setattr__('starttime', dt0)
            proc.__setattr__('endtime', dt1)
            stt = prepare_data(st, proc)
            run_process(stt, proc)
    else:
        stt = prepare_data(st, proc)
        run_process(stt, proc)

    return

# Helper functions

def array_processing(st, proc, verbose=False):
    stationtable = export_stationtable(proc)
    data_files = export_datafiles(st)
    data_files = " ".join(data_files)

    if proc.method == 'ccts':
        print (' - Cross-correlation Trace Stacking (CCTS) analysis...')
        oversampling = 1
        cmd  = f'{proc.method} {stationtable} ' \
        f'{proc.binsize} {proc.overlap} {oversampling} ' \
        f'{proc.theta[0]} {proc.theta[1]} {proc.theta[2]} ' \
        f'{proc.c_trace[0]} {proc.c_trace[1]} {proc.c_trace[2]} ' \
        f'{proc.grid_type} {proc.grid_return} sac {data_files}'

    elif proc.method == 'freqfisher':
        print (' - FK Fisher analysis...')
        cmd  = f'{proc.method} {stationtable} ' \
               f'{proc.binsize} {proc.overlap} ' \
               f'{proc.freq[0]} {proc.freq[1]} ' \
               f'{proc.freq[2]} {proc.freq[3]} ' \
               f'{proc.theta[0]} {proc.theta[1]} {proc.theta[2]} ' \
               f'{proc.c_trace[0]} {proc.c_trace[1]} {proc.c_trace[2]} ' \
               f'{proc.grid_type} sac {data_files}'

    elif proc.method == 'tdoa':
        print (' - Time Delay of Arrival (TDOA) analysis...')
        oversampling = 1
        cmd  = f'{proc.method} {stationtable} ' \
        f'{proc.binsize} {proc.overlap} ' \
        f'{oversampling} sac {data_files}'

    elif proc.method == 'tfreqfisher':
        print (' - FK+ Fisher analysis...')
        td_file = '{}_{}.dat'.format(proc.timedomain_detector_method,
                                     proc.starttime.strftime('%Y%m%d'))
        shutil.copyfile(td_file, 'timefisher.dat')
        cmd  = f'{proc.method} {stationtable} ' \
               f'{proc.binsize} {proc.overlap} ' \
               f'{proc.freq[0]} {proc.freq[1]} ' \
               f'{proc.freq[2]} {proc.freq[3]} ' \
               f'sac {data_files}'

    elif proc.method == 'timefisher':
        print (' - Time-domain Fisher analysis...')
        cmd  = f'{proc.method} {stationtable} ' \
        f'{proc.binsize} {proc.overlap} ' \
        f'{proc.theta[0]} {proc.theta[1]} {proc.theta[2]} ' \
        f'{proc.c_trace[0]} {proc.c_trace[1]} {proc.c_trace[2]} ' \
        f'{proc.grid_type} {proc.grid_return} sac {data_files}'
    
    if verbose:
        print(cmd)
    os.system(cmd)

    processor_output = f'{proc.method}.dat'
    timestamp = proc.starttime.strftime('%Y%m%d')
    results_file = '{}_{}.dat'.format(proc.method,
                                      timestamp)
    print (' - Writing output to [ %s ]' % results_file)
    os.rename(processor_output, results_file)
    try:
        bestbeam_file = 'bestbeam_{}.sac'.format(timestamp)
        os.rename('bestbeam.sac', bestbeam_file)
    except:
        pass
    return

def run_process(st, proc):
    """
    Wrapper function to call processes. Only when data is present
    """
    if len(st) > 0:
        print('')
        print(st)
        print('')
        array_processing(st, proc)
    return

def export_datafiles(st):
    """
    Function to write files for array processing as SAC files
    """
    data_files = []

    for tr in st:
        fid = '{}.SAC'.format(tr.id)
        tr.write(fid, 'SAC')
        data_files.append(fid)
    return data_files

def export_stationtable(proc):
    """
    Wrapper function to export inventory to stationtable
    """
    df = inv_to_df(proc.inventory, compute_offset=True)
    fid_stationtable = proc.name+'.dat'
    df_to_ascii(df, fid_stationtable, print_offset=True)
    return fid_stationtable

def prepare_data(st, proc):
    print ('Processing timespan {} - {}'.format(proc.starttime, proc.endtime))

    stt = trim_stream(st, proc)
    if len(stt) > 0:
        # Throwing away bad traces
        for tr in stt:
            tr.data = tr.data.astype(np.float)
            if tr.stats.npts == 1:
                stt.remove(tr)
    
        print (' - Merging data...')
        stt.merge(fill_value='interpolate')
        stt.sort(keys=['station','starttime'])

        print (' - Detrending data...')
        stt.detrend()

        fmin = proc.freq[0]
        fmax = proc.freq[1]
        print (' - Filtering data [{:.2f} - {:.2f}] Hz'.format(fmin, fmax))
        if proc.deconvolution == 'response':
            print ('  - Deconvolving instrument response...')
            stt.attach_response(proc.inventory)
            pre_filt = [0.9*fmin, fmin, fmax, 1.1*fmax]
            stt.remove_response(pre_filt=pre_filt)

        elif proc.deconvolution == 'gain':
            print ('  - Correcting for instrument gain...')
            for tr in stt:
                gain = get_gain(proc.inventory, tr,
                                starttime=proc.starttime,
                                endtime=proc.endtime)
                tr.data *= 1.0/gain

            stt.filter('bandpass', corners=4, zerophase=True,
                    freqmin=fmin, freqmax=fmax)

        if proc.sampling_rate:
            print (' - Anti-alias filter and resampling data...')
            for tr in stt:
                fs = tr.stats.sampling_rate
                if  fs != proc.sampling_rate:
                    tr.interpolate(proc.sampling_rate)
    return stt

def trim_stream(st, proc):
    """
    Function to trim dataset to start end endtime
    """
    stt = Stream()
    for tr in st:
        dt = tr.stats.delta
        stt += tr.copy().trim(proc.starttime, proc.endtime-dt)
    return stt

def comma_list(s):
    return s.split(',')

def parameter_list(s):
    return list(map(float, s.split('/')))

def process_cli(argv):
    parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument('-n', '--name',
        type=str, metavar='Array name'
        )
    parser.add_argument('-m', '--method',
        type=str, choices=bf_methods
        )
    parser.add_argument('-inv', '--inventory',
        type=str, metavar='StationXML inventory (.xml)'
        )
    parser.add_argument('-st', '--starttime',
        type=lambda t: UTCDateTime(t),
        metavar='Start time (ObsPy compatible string)'
        )
    parser.add_argument('-et', '--endtime',
        type=lambda t: UTCDateTime(t),
        metavar='End time (ObsPy compatible string)'
        )
    parser.add_argument('-bs', '--binsize',
        type=float, metavar='Time-bin (s)'
        )
    parser.add_argument('-ov', '--overlap',
        type=int, metavar='Overlap window (%)'
        )
    parser.add_argument('-fr', '--freq',
        type=parameter_list, metavar='fmin/fmax (td) or fmin/fmax/fstep/fav (fd)'
        )
    parser.add_argument('-dr', '--data_root',
        type=str, metavar='Root directory SDS structure'
        )

    # Optional arguments
    parser.add_argument('-th', '--theta',
        nargs='?', type=parameter_list, metavar='th_min/th_max/dth'
        )
    parser.add_argument('-ct', '--c_trace',
        nargs='?', type=parameter_list, metavar='c_min/c_max/dc'
        )
    parser.add_argument('-gt', '--grid_type',
        nargs='?', type=str, choices=['tele','local']
        )
    parser.add_argument('-gr', '--grid_return',
        nargs='?', type=str, choices=['max','all']
        )
    parser.add_argument('-el','--skip_elements',
        nargs='?', type=comma_list,
        metavar='Stations to skip: <NET.STA.LOC.CHAN>, <NET.STA.LOC.CHAN>'
        )
    parser.add_argument('-decon', '--deconvolution',
        nargs='?', type=str, choices=['gain','response']
        )
    parser.add_argument('-fs', '--sampling_rate',
        nargs='?', type=float,
        metavar='Sample rate of processing (Hz)'
        )
    parser.add_argument('-tdm', '--timedomain_detector_method',
        nargs='?', type=str, choices=['ccts', 'tdoa', 'timefisher']
        )
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    return args

if __name__ == "__main__":
   main(sys.argv[1:])
