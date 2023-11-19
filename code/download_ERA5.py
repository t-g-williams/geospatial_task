'''
Download data from the ERA5 portal
'''
import cdsapi
import pandas as pd
import numpy as np

base_dir = 'D:/My Drive/job_applications/2023_industry/climate_farmers/geospatial_task/'

# format some data for the query
# bbox: read from file
bbox = pd.read_csv(f'{base_dir}data/region_extent.csv')
bbox_list = [bbox.ymax[0], bbox.xmin[0], bbox.ymin[0], bbox.xmax[0]]
# years: 2000 to 2022
years = list(np.arange(2000,2023).astype(str))

# run the query to download the data
c = cdsapi.Client()
c.retrieve(
    'reanalysis-era5-land-monthly-means',
    {
        'product_type': 'monthly_averaged_reanalysis',
        'variable': [
            '2m_temperature', 'potential_evaporation', 'soil_temperature_level_1',
            'total_evaporation', 'total_precipitation',
            ],
        'year': years,
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',
        'area': bbox_list,
        'format': 'netcdf',
    },
    f'{base_dir}data/biophysical/ERA5_download.nc'
    )