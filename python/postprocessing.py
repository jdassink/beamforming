"""
Beamforming post processing helper functions

.. module:: postprocessing

:author:
    Jelle Assink (jelle.assink@knmi.nl)

:copyright:
    2021, Jelle Assink

:license:
    This code is distributed under the terms of the
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.en.html)
"""
import xarray as xr
import pandas as pd
import numpy as np

import pyproj
from datetime import datetime
from scipy.stats import binned_statistic_2d

def bin_xyz_values(ds, config, dx, dy, statistic='mean'):
    """
    Function to grid x, y, z data from three columns using SciPy.
    Note that time data must be converted to a numerical value first.

    Parameters
    ----------
    ds : `XArray Dataset`
        Dataset containing x, y, z input data

    config : `dict`
        dictionary mapping x, y, z columns to Dataset variables
        e.g. config=dict(x='time', y='baz', z='snr')

    dx, dy : `float`
        x, y binsize

    statistic : `mean`, 
        Statistic that is used as bin value. See `scipy.stats.binned_statistic_2d`

    Returns
    -------
    ds : `xarray Dataset`
        Dataset with binned data
    """
    if config['x'] == 'time':
        x = ds[config['x']].values.astype(np.timedelta64) / np.timedelta64(1, 's')
    else:
        x = ds[config['x']].values
    y = ds[config['y']].values
    z = ds[config['z']].values

    x_bins = np.arange(np.nanmin(x), np.nanmax(x)+dx, dx)
    y_bins = np.arange(np.nanmin(y), np.nanmax(y)+dy, dy)

    (Z, X, Y, _) = binned_statistic_2d(x, y, z, statistic=statistic,
                                       bins=[x_bins, y_bins])
    if config['x'] == 'time':
        X = [ datetime.utcfromtimestamp(ts) for ts in X ]
    
    x_coord = '{}_bin'.format(config['x'])
    y_coord = '{}_bin'.format(config['y'])

    da = xr.DataArray(
        data = Z.T,
        dims=[y_coord, x_coord],
        coords={ x_coord : ([ x_coord ], X[:-1]),
                 y_coord : ([ y_coord ], Y[:-1])} )
    ds = da.to_dataset(name=config['y'])
    name = '{} {}'.format(config['y'], statistic)
    ds.attrs = dict(long_name=name,
                    standard_name=name,
                    standard_description=statistic)
    return ds


def compute_weighted_time_average(ds, minimum_count=1):
    """
    Compute weighted average for binned data

    Parameters
    ----------
    ds : `XArray Dataset`
        Dataset containing binned input data

    minimum_count: `int`
        Minimum number of detections in bin to be considered for average

    Returns
    -------
    ds : `xarray Dataset`
        Dataset with binned data
    """
    dsc = ds.copy(deep=True)

    vars = list(ds.data_vars)
    y_coord = list(ds.coords)[1]

    dswa = []
    # compute weighted average
    for var in vars:
        ds_ = dsc[var].where(ds[var] >= minimum_count, other=0)
        da = np.matmul(ds_.values.T, ds_[y_coord].values)
        sumval = ds_.sum(dim=y_coord, skipna=True)
        sumval = sumval.where(sumval > 0)
        da /= sumval
        varname = f'{var}_wgt_av'
        dswa.append(da.to_dataset(name=varname))

    return xr.merge(dswa)

def freqfisher_to_dataset(ff_file, t0):
    """
    Function to read time-domain Fisher analysis results and convert to XArray Dataset

    Parameters
    ----------
    ff_file : `str`
        File name of ASCII file containing time-domain Fisher analysis results

    t0 : `datetime`
        Datetime object with reference time of processing file

    Returns
    -------
    ds : `dictionary`
        Dictionary of `XArray.Dataset` objects
    """
    colnames = ['time','freq', 'fratio', 'baz', 'app_vel', 'n_instr', 'psd']

    df = pd.read_csv(ff_file,
                     names=colnames,
                     delimiter='\s+',
                     engine='c')

    df['time'] = pd.to_datetime(df.time, unit='s', origin=t0)
    
    n_freq = df.freq.nunique()
    n_bins = df.time.nunique()
    freq = df.freq.unique()
    time = df.time.unique()

    # rework values
    snr = (df.fratio.values - 1) / df.n_instr.values
    snr = snr.reshape(n_bins, n_freq)
    snr = np.sqrt(snr)
    fratio = df.fratio.values.reshape(n_bins, n_freq)
    app_vel = df.app_vel.values.reshape(n_bins, n_freq)
    baz = df.baz.values.reshape(n_bins, n_freq)
    psd = df.psd.values.reshape(n_bins, n_freq)

    reference = ('Smart, E. and Flinn, E. A. (1971), '
                'Fast Frequency‐Wavenumber Analysis and Fisher Signal Detection '
                'in Real‐Time Infrasonic Array Data Processing. '
                'Geophysical Journal of the Royal Astronomical Society, 26: 279-284. '
                'http://dx.doi.org/10.1111/j.1365-246X.1971.tb03401.x')

    general_attrs = dict(references=reference,
                        software='https://github.com/jdassink/beamforming',
                        institution='KNMI',
                        author='Jelle Assink',
                        email='assink@knmi.nl')

    snr_attrs = dict(long_name='SNR', units='-')
    fratio_attrs = dict(long_name='Fisher ratio', units='-')
    baz_attrs = dict(long_name='Back azimuth', units='degrees')
    vel_attrs = dict(long_name='Apparent velocity', units='m/s')
    psd_attrs = dict(long_name='Power spectral density', units='Pa^2/Hz')

    dims = ["freq", "time"]

    ds = xr.Dataset(data_vars=dict(snr=(dims, snr.T, snr_attrs),
                                fratio=(dims, fratio.T, fratio_attrs),
                                app_vel=(dims, app_vel.T, vel_attrs),
                                baz=(dims, baz.T, baz_attrs),
                                psd=(dims, psd.T, psd_attrs)),
                    coords=dict(freq=(["freq"], freq),
                                time=(["time"], time)),
                    attrs=general_attrs)

    return ds

def timefisher_to_dataset(tf_file, t0):
    """
    Function to read time-domain Fisher analysis results and convert to XArray Dataset

    Parameters
    ----------
    tf_file : `str`
        File name of ASCII file containing time-domain Fisher analysis results

    t0 : `datetime`
        Datetime object with reference time of processing file

    Returns
    -------
    ds : `dictionary`
        Dictionary of `XArray.Dataset` objects
    """
    colnames = ['time', 'fratio',
                'baz', 'app_vel',
                'sx', 'sy',
                'prms', 'p2p', 'n_instr',
                'f_center', 'f_bw', 'f_rms']

    df = pd.read_csv(tf_file,
                     names=colnames,
                     usecols=[0, 1, 2, 3, 6, 7, 8, 9],
                     delimiter='\s+',
                     engine='c')

    df['time'] = pd.to_datetime(df.time, unit='s', origin=t0)
    df['snr'] = np.sqrt((df['fratio'] - 1)/df['n_instr'])
    df = df.set_index('time')
    return df.to_xarray()

def tdoa_to_dataset(tf_file, t0):
    """
    Function to read time-delay-of-arrival analysis results and convert to XArray Dataset

    Parameters
    ----------
    tf_file : `str`
        File name of ASCII file containing time-delay-of-arrival analysis results

    t0 : `datetime`
        Datetime object with reference time of processing file

    Returns
    -------
    ds : `dictionary`
        Dictionary of `XArray.Dataset` objects
    """
    colnames = ['time', 'fratio',
                'baz', 'app_vel',
                'sx', 'sy',
                'prms', 'p2p', 'n_instr',
                'f_center', '-',
                'mccm', 'lsq_error', 'sigma_tau',
                'q_value',
                'sigma_baz', 'sigma_vel']

    df = pd.read_csv(tf_file,
                     names=colnames,
                     usecols=[0, 1, 2, 3, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16],
                     delimiter='\s+',
                     engine='c')

    df['time'] = pd.to_datetime(df.time, unit='s', origin=t0)
    df['snr'] = np.sqrt((df['fratio'] - 1)/df['n_instr'])
    df = df.set_index('time')
    return df.to_xarray()

def inverse_transform(start, end):
    """
    Helper function. Returns (back)azimuth and distance on WGS-84 ellipsoid

    Parameters
    ----------
    start : `dict`
        Dictionary w/ latitude ('lat') and longitude ('lon') of start point

    end : `dict`
        Dictionary w/ latitude ('lat') and longitude ('lon') of end point

    Returns
    -------
    azi : `float`
        Azimuth (between 0-360 deg) from start to end point

    bazi : `float`
        Back azimuth (between 0-360 deg) from end to start point

    dist : `float`
        Distance (in meters) between start end end point
    """
    g = pyproj.Geod(ellps='WGS84')
    (azi, bazi, dist) = g.inv(start['lon'], start['lat'],
                                end['lon'], end['lat'], radians=False)
    return(azi%360, bazi%360, dist)
