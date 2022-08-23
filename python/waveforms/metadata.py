"""
Metadata interface functions

.. module:: metadata

:author:
    Jelle Assink (jelle.assink@knmi.nl)

:copyright:
    2021, Jelle Assink

:license:
    This code is distributed under the terms of the
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.en.html)
"""
from obspy import read_inventory
from obspy.clients.fdsn import Client as fdsn_client
from obspy.clients.fdsn import RoutingClient
from obspy.core.inventory.inventory import Inventory
from obspy.signal.util import util_geo_km
from obspy.geodetics.base import gps2dist_azimuth

from pandas import DataFrame
from lxml import etree as ET
import numpy as np


class param_object(object):
    pass


def compute_aperture(coords):
    aperture = 0

    n_sta = len(coords)
    for i in range(0, n_sta):
        for j in range(i+1, n_sta):
            (d, _, _) = gps2dist_azimuth(coords[i][0], coords[i][1],
                                         coords[j][0], coords[j][1])
            if d > aperture:
                aperture = d
    return aperture


def compute_offsets(cha, ref):
    (x, y) = util_geo_km(ref.longitude, ref.latitude,
                         cha.longitude, cha.latitude)
    return (x, y)


def df_to_ascii(df, fid_out, print_offset=False, verbose=False):
    """
    Write out stationtable as formatted ASCII table
    """
    if verbose:
        print(f'Writing to {fid_out}')

    with open(fid_out, 'w') as f:
        for _, r in df.iterrows():
            seed_id = f'{r.network}.{r.station}.{r.location}.{r.channel}'
            line = f'{seed_id:>18s} '

            if print_offset:
                line += f'{r.deast:17.4f} {r.dnorth:17.4f} '
                line += f'{r.elevation:17.4f} {r.edepth:17.4f} '
                line += f'{r.latitude:9.6f} {r.longitude:9.6f} '
                line += f'{r.starttime:8d} {r.endtime:8d} '
                line += f'{r.gain:15.4e}\n'
            f.write(line)
    return


def get_array_reference(inv, reference='first'):
    ref = param_object()

    if reference == 'first':
        ref.latitude = inv[0][0][0].latitude
        ref.longitude = inv[0][0][0].longitude

    elif reference == 'center':
        data = []
        for seed_id in inv.get_contents()['channels']:
            data.append(inv.get_coordinates(seed_id))
        data = DataFrame(data).drop_duplicates()
        ref.latitude = data['latitude'].mean()
        ref.longitude = data['longitude'].mean()
    return ref


def get_gain(inv, tr, starttime=None, endtime=None):
    invs = inv.select(network=tr.stats.network,
                      station=tr.stats.station,
                      location=tr.stats.location,
                      channel=tr.stats.channel,
                      starttime=starttime,
                      endtime=endtime)
    return invs[0][0][0].response.instrument_sensitivity.value


def get_metadata(data_source, network, station, location, channel,
                 starttime, endtime, **kwargs):
    """
    source type options:

    local - local stationxml file
    fdsn - typical dataselect
    routing - 'iris-federator' or 'eida-routing'
    """
    invselect = locals()
    del(invselect['data_source'])
    del(invselect['kwargs'])
    [(key, value)] = data_source.items()

    if key == 'local':
        inv = read_inventory(value, 'STATIONXML')
        inv = inv.select(**invselect)

    elif key == 'fdsn':
        cl = fdsn_client(value)
        inv = cl.get_stations(**invselect, **kwargs, level='response')

    elif key == 'routing':
        cl = RoutingClient(value)
        inv = cl.get_stations(**invselect, **kwargs, level='response')

    else:
        print('meta-data source should be "local" or "fdsn"')
        inv = []
    return inv


def inv_to_coords(inv):
    """
    Convert inventory with items to Numpy array for array
    response calculations
    """
    n_stations = len(inv.get_contents()['stations'])
    coords = np.zeros((n_stations, 3))
    i = 0

    for i in range(0, n_stations):
        coords[i, 0] = inv[0][i].latitude
        coords[i, 1] = inv[0][i].longitude
        coords[i, 2] = inv[0][i].elevation

    return coords


def inv_to_df(inv, compute_offset=False, offset_reference='first'):
    data = []

    for net in inv:
        for sta in net:
            for cha in sta:
                t0 = parse_utcdatetime(cha.start_date)
                t1 = parse_utcdatetime(cha.end_date)
                try:
                    gain = cha.response.instrument_sensitivity.value
                except ValueError as e:
                    print(e)
                    gain = -1
                item = dict(network=net.code,
                            station=sta.code,
                            location=cha.location_code,
                            channel=cha.code,
                            latitude=cha.latitude,
                            longitude=cha.longitude,
                            elevation=cha.elevation,
                            edepth=cha.depth,
                            starttime=t0,
                            endtime=t1,
                            gain=gain)

                if compute_offset:
                    ref = get_array_reference(inv, offset_reference)
                    (x, y) = compute_offsets(cha, ref)
                    item['deast'] = x*1e3
                    item['dnorth'] = y*1e3

                data.append(item)
    df = DataFrame(data).drop_duplicates()
    return df


def inv_to_dtk_xml(inv, array_name):
    """
    Function to convert StationXML to XML tree object
    that is compatible with DTK-GPMCC
    """
    ns = {None: 'http://www.fdsn.org/xml/station/1',
          'dtk': 'http://www.fdsn.org/xml/dtk/1'}

    date_fmt = '%Y-%m-%dT%H:%M:%S'
    root = ET.Element('FDSNStationXML', nsmap=ns, schemaVersion=u"1")

    for net in inv:
        network = ET.SubElement(root, 'Network', code=f"{net.code}")
        nsta = ET.SubElement(network, 'TotalNumberStations')
        nsta.text = f'{len(inv.get_contents()["networks"])}'

        ref_sta = net[0]
        t0 = ref_sta.start_date.strftime(date_fmt)
        try:
            t1 = ref_sta.end_date.strftime(date_fmt)
        except ValueError as e:
            print(e)
            t1 = None
        station = ET.SubElement(network, 'Station',
                                code=f'{array_name}',
                                startDate=f'{t0}',
                                endDate=f'{t1}')

        latitude = ET.SubElement(station, 'Latitude')
        latitude.text = f'{ref_sta.latitude}'
        longitude = ET.SubElement(station, 'Longitude')
        longitude.text = f'{ref_sta.longitude}'
        elevation = ET.SubElement(station, 'Elevation')
        elevation.text = f'{ref_sta.elevation}'

        ncha = ET.SubElement(station, 'TotalNumberChannels')
        ncha.text = f'{len(inv.get_contents()["channels"])}'

        for item in inv.get_contents()['channels']:
            (net, sta, loc, cha) = item.split('.')
            invs = inv.select(network=net, station=sta,
                              location=loc, channel=cha)

            for _cha in invs[0][0]:
                t0 = _cha.start_date.strftime(date_fmt)
                try:
                    t1 = _cha.end_date.strftime(date_fmt)
                except ValueError as e:
                    print(e)
                    t1 = None
                channel = ET.SubElement(station, 'Channel',
                                        code=f'{cha}',
                                        name=f'{sta}',
                                        locationCode=f'{loc}',
                                        startDate=f'{t0}',
                                        endDate=f'{t1}')

                latitude = ET.SubElement(channel, 'Latitude')
                latitude.text = f'{_cha.latitude}'
                longitude = ET.SubElement(channel, 'Longitude')
                longitude.text = f'{_cha.longitude}'
                elevation = ET.SubElement(channel, 'Elevation')
                elevation.text = f'{_cha.elevation/1e3}'
                depth = ET.SubElement(channel, 'Depth')
                depth.text = f'{_cha.depth}'
                srate = ET.SubElement(channel, 'SampleRate')
                srate.text = f'{_cha.sample_rate}'

                resp = ET.SubElement(channel, 'Response')
                _resp = _cha.response
                sens = ET.SubElement(resp, 'InstrumentSensitivity')
                _sens = _resp.instrument_sensitivity
                out_units = ET.SubElement(sens, 'OutputUnits')
                out_units.text = _sens.output_units
                sens_value = ET.SubElement(sens, 'Value')
                sens_value.text = f'{_sens.value}'
                sens_freq = ET.SubElement(sens, 'Frequency')
                sens_freq.text = f'{_sens.frequency}'
                in_units = ET.SubElement(sens, 'InputUnits')
                in_units_name = ET.SubElement(in_units, 'Name')
                in_units_name.text = _sens.input_units

    tree = ET.ElementTree(root)
    return tree


def parse_utcdatetime(time):
    """
    Convenience function to format UTCDateTime in year-julday format
    """
    try:
        val = int(time.strftime('%Y%j'))
    except ValueError as e:
        print(e)
        val = -1
    return val


def select_inventory(stream, inv, starttime=None, endtime=None):
    """
    Convenience function to select inventory from active traces only
    """
    invs = Inventory()
    for tr in stream:
        invs += inv.select(network=tr.stats.network,
                           station=tr.stats.station,
                           location=tr.stats.location,
                           channel=tr.stats.channel,
                           starttime=starttime,
                           endtime=endtime)
    return invs
