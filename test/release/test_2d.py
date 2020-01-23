import datetime
import hm.api as hm
import xarray
import pytest

# os.chdir("test/release")

def test_set_domain():
    modeltime = hm.set_modeltime(
        datetime.datetime(2003,6,1),
        datetime.datetime(2003,6,30),
        datetime.timedelta(days=1)
    )
    hm.set_domain(
        'test_data/ghana_landmask_0pt25degree.tif',
        modeltime
    )
    
def test_open_hmdataarray():
    modeltime = hm.set_modeltime(
        datetime.datetime(2003,6,1,0,0),
        datetime.datetime(2003,6,10,0,0),
        datetime.timedelta(days=1)
    )
    domain = hm.set_domain(
        'test_data/ghana_landmask_0pt25degree.tif',
        modeltime
    )
    da = hm.open_hmdataarray(
        'test_data/AgMERRA_2003_tavg.nc4',
        'tavg',
        domain
    )
    da.load()
    
    mask = da._domain.mask
    da = da._data.transpose()
    da = da.sel(time=pd.Timestamp(2003,6,1), method='nearest')
    
import xarray as xr
mask = xr.open_rasterio('test_data/ghana_landmask_0pt25degree.tif')
mask = mask.sel(band=1)
da = xr.open_dataset('test_data/AgMERRA_2003_tavg.nc4')['tavg']
# da2 = da.transpose('latitude','longitude','time')
da = xr.load_dataset('/data/WFDEI/PSurf_daily_WFDEI/PSurf_daily_WFDEI_200001.nc', decode_times=False)


import netCDF4 as nc
x = nc.Dataset('test_data/AgMERRA_2003_tavg.nc4')
