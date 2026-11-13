from obspy.core.inventory import (
    Inventory,
    util,
    Network,
    Station,
    Channel,
    Site,
    Response,
    Comment
)

import obspy
from obspy.core.inventory.util import Frequency
import argparse
from pathlib import Path
import pandas as pd

'''
Create stationxml files for COSZO stations. Separate files are created for each location code. E.g., OO_CZMID_30.xml, OO_CZMID_31.xml

Requirements: COSZO Instrument Data Channels.xlsx

Example usage:

python create_coszo_metadata.py --metadata "../param/COSZO Instrument Data Channels.xlsx" --output ../output/xml --station-name AXCC1 --location-code 40
(the spreadsheet is the COSZO team's channel table; it is not part of this public repository)
'''

parser = argparse.ArgumentParser(description='Create StationXML from COSZO metadata spreadsheet.')
parser.add_argument('--metadata', required=True, help='Path to the metadata spreadsheet (.xlsx)')
parser.add_argument('--station-name', required=True, help='Station code to filter on, e.g. CZMID')
parser.add_argument('--location-code', required=True, help='Location code to filter on, e.g. 30')
parser.add_argument(
    '--exclude-channels',
    nargs='*',
    default=[],
    help='Channel codes to leave out, e.g. channels without data yet'
)
parser.add_argument(
    '--output',
    required=False,
    default='.',
    help='Output XML file path or output directory. If a directory is given, the filename is generated automatically.'
)
args = parser.parse_args()

metadata = args.metadata
station_name = args.station_name
location_code = str(args.location_code)
output_path = Path(args.output)

# Returns None if field in metadata spreadsheet is left blank
def blank_to_none(x):
    return None if pd.isna(x) else x

# Dataframe with network information
df_n = pd.read_excel(metadata, skiprows=1, header=0, sheet_name='Network Metadata')

# Dataframe with all channels
df = pd.read_excel(metadata, skiprows=[0,1,2,3,5], header=0, sheet_name='Channel Metadata')

# Filter dataframe by station name and location code
df = df[df['Station Code'] == station_name]
df = df[pd.to_numeric(df['Channel Location Code'], errors='coerce') == pd.to_numeric(location_code)]
if args.exclude_channels:
    df = df[~df['Channel Code'].isin(args.exclude_channels)]
    print('Excluded channels:', ', '.join(args.exclude_channels))
print('Number channels in dataframe =', len(df))

if df.empty:
    raise ValueError(f'No channels found for station {station_name!r} and location code {location_code!r}')

inv = Inventory(
    networks=[],
    source=df_n['Source'].iloc[0],
) 

if pd.isna(df_n['Network End Date'].iloc[0]) == True:
    end_date=None
else:
    end_date = obspy.UTCDateTime(df_n['Network End Date'].iloc[0])

# Get network information
net = Network(
    code = df_n['Network Code'].iloc[0],
    stations = [],
    description = df_n['Network Description'].iloc[0],
    start_date = obspy.UTCDateTime(df_n['Network Start Date'].iloc[0]),
    end_date = end_date,
    restricted_status=df_n['Network Status'].iloc[0],
)

# Get station information
if pd.isna(df['Station End Date'].iloc[0]) == True:
    end_date = None
else:
    end_date = obspy.UTCDateTime(df['Station End Date'].iloc[0])

sta = Station(
    code = df['Station Code'].iloc[0],
    start_date = obspy.UTCDateTime(df['Station Start Date'].iloc[0]),
    end_date = end_date,
    latitude = df['Station Latitude'].iloc[0],
    longitude = df['Station Longitude'].iloc[0],
    elevation = df['Station Elevation'].iloc[0],
    site=Site(name = df['Station Site'].iloc[0]),
    restricted_status = df['Station Status'].iloc[0],
)

# Get channel information 
cha = []
channels = []

for i in range(len(df)):

    sample_rate = float(df['Channel Sample Rate (Hz)'].iloc[i])

    # Get instrument sensitivity
    sensitivity = obspy.core.inventory.response.InstrumentSensitivity(
        value = float(df['Response Instrument Sensitivity Value'].iloc[i]),
        frequency = df['Response Instrument Sensitivity Frequency'].iloc[i],
        input_units = df['Response Instrument Sensitivity Input Units'].iloc[i],
        input_units_description =df['Response Input Units Description'].iloc[i],
        output_units = df['Response Instrument Sensitivity Output Units'].iloc[i],
        output_units_description = df['Response Output Units Description'].iloc[i],
    )

    # Get instrument response
    response_stage = obspy.core.inventory.response.PolesZerosResponseStage(
        stage_sequence_number = 1,
        stage_gain = float(df['Response Instrument Sensitivity Value'].iloc[i]),
        stage_gain_frequency = df['Response Instrument Sensitivity Frequency'].iloc[i],
        input_units = df['Response Instrument Sensitivity Input Units'].iloc[i],
        input_units_description = df['Response Input Units Description'].iloc[i],
        output_units = df['Response Instrument Sensitivity Output Units'].iloc[i],
        output_units_description = df['Response Output Units Description'].iloc[i],
        pz_transfer_function_type="DIGITAL (Z-TRANSFORM)",

        normalization_frequency=1.0,
        normalization_factor=1.0,
        zeros=[],
        poles=[],
        resource_id=None,
        resource_id2=None,
        name="Stage01",
        decimation_input_sample_rate=Frequency(sample_rate),
        decimation_factor=1,
        decimation_offset=0,
        decimation_delay=Frequency(0.0),
        decimation_correction=Frequency(0.0),
    )

    res = Response(
        instrument_sensitivity=sensitivity, response_stages=[response_stage]
    )

    if pd.isna(df['Channel End Date'].iloc[0]) == True:
        end_date=None
    else:
        end_date = blank_to_none(obspy.UTCDateTime(df['Channel End Date'].iloc[i]))

    cha = Channel(
        code = df['Channel Code'].iloc[i],
        location_code = str(df['Channel Location Code'].iloc[i]),
        start_date=obspy.UTCDateTime(df['Channel Start Date'].iloc[i]),
        end_date=end_date,
        latitude = df['Channel Latitude'].iloc[i],
        longitude = df['Channel Longitude'].iloc[i],
        elevation = df['Channel Elevation'].iloc[i],
        depth = df['Channel Depth'].iloc[i],
        azimuth = blank_to_none(df['Channel Azimuth'].iloc[i]),
        dip = blank_to_none(df['Channel Dip'].iloc[i]),
        sample_rate = (df['Channel Sample Rate (Hz)'].iloc[i]),
        description = df['Channel Description'].iloc[i],
        sensor=util.Equipment(description=df['Channel Sensor Description'].iloc[i]),
        response=res,
    )

    # Blank comment cells come through as NaN; skip them so no "nan" comment is written.
    comment_value = blank_to_none(df['Channel Comment Value'].iloc[i])
    if comment_value is not None:
        cha.comments.append(Comment(value=str(comment_value)))

    channels.append(cha)

# Tie it all together.
inv.networks.append(net)

# Added this check to avoid TypeError None <= datetime(...). MT
def leq(a, b):
    if a is None and b is None: 
        return True
    if a is None:               
        return False
    if b is None:              
        return True
    return a <= b

# Added this check to avoid TypeError None <= datetime(...). MT
def geq(a, b):
    if a is None and b is None: 
        return True
    if a is None:               # a is −∞? If you want −∞ for start, return False
        return False
    if b is None:               # b is −∞ → a >= −∞ is always True
        return True
    return a >= b

if leq(sta.end_date, net.end_date) and geq(sta.start_date, net.start_date):
    net.stations.append(sta)
    for cha in channels:
        if leq(cha.end_date, sta.end_date) and geq(cha.start_date, sta.start_date):
            sta.channels.append(cha)

for s in net.stations:
    print("Station:", s.code, s.start_date, s.end_date, "channels:", len(s.channels))

default_filename = f"{df_n['Network Code'].iloc[0]}_{df['Station Code'].iloc[0]}_{location_code}.xml"
if output_path.suffix.lower() == '.xml':
    output_file = output_path
else:
    output_file = output_path / default_filename

output_file.parent.mkdir(parents=True, exist_ok=True)

inv.write(
    str(output_file),
    format="STATIONXML",
    validate=True,
)
print(f'Wrote {output_file}')
