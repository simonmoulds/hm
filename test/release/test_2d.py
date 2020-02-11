import datetime
import hm.api as hm
import xarray
import pytest
import pandas as pd

# os.chdir("test/release")
# modeltime = hm.set_modeltime(
#     datetime.datetime(2003,6,1,0,0),
#     datetime.datetime(2003,6,10,0,0),
#     datetime.timedelta(days=1)
# )
# domain = hm.set_domain(
#     'test_data/ghana_landmask_0pt25degree.tif',
#     modeltime
# )
# da = hm.open_hmdataarray(
#     'test_data/AgMERRA_2003_tavg.nc4',
#     'tavg',
#     domain
# )
# x = da.select(time=pd.Timestamp(2003,6,4), method='nearest')
# x = da.select(time=slice(pd.Timestamp(2003,6,4), pd.Timestamp(2003,6,6)))

# def match(x, table):
#     table_sorted = np.argsort(table)
#     x_pos = np.searchsorted(table[table_sorted], x)
#     return table_sorted[x_pos]

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
    da.select(time=pd.Timestamp(2003,6,1), method='nearest')
    # da.load()    
    # mask = da._domain.mask
    # da = da._data.transpose()
    # da = da.sel(time=pd.Timestamp(2003,6,1), method='nearest')
    
# import xarray as xr
# mask = xr.open_rasterio('test_data/ghana_landmask_0pt25degree.tif')
# mask = mask.sel(band=1)
# da = xr.open_dataset('test_data/AgMERRA_2003_tavg.nc4')
# # da2 = da.transpose('latitude','longitude','time')
# da = xr.load_dataset('/data/WFDEI/PSurf_daily_WFDEI/PSurf_daily_WFDEI_200001.nc', decode_times=False)


# import netCDF4 as nc
# import xarray as xr
# x = nc.Dataset('test_data/AgMERRA_2003_tavg.nc4', mode='r')
# ds = xr.open_dataset('test_data/AgMERRA_2003_tavg.nc4')
# temp = ds.variables['tavg'].values
# x = ds.variables['longitude'].values
# y = ds.variables['latitude'].values
# tm = ds.variables['time'].values
# # dsn = xr.Dataset({'tavg' : (['time','latitude','longitude'], temp)},
# #                  coords={'longitude' : x,
# #                          'latitude' : y,
# #                          'time' : tm})
# dsn = xr.Dataset({'tavg' : ['time','latitude','longitude']},
#                  coords={'longitude' : x,
#                          'latitude' : y,
#                          'time' : tm})


# modeltime = hm.set_modeltime(
#     datetime.datetime(2003,6,1,0,0),
#     datetime.datetime(2004,6,10,0,0),
#     datetime.timedelta(days=1)
# )

# coords = {}
# with xr.open_dataset('test_data/AgMERRA_2003_tavg.nc4') as ds:
#     for dim in x.dimensions:
#         coords[dim] = ds[dim].values


# a = np.array([1,2,3,4,5,6,7,8,9,10])
# b = np.array([1,5,1,10,10,6,4,4])
# ab,a_ind,b_ind = np.intersect1d(a,b,return_indices=True)

# a = np.random.randint(0,100000,(100000,))
# b = np.random.randint(0,100000,(100000,))
# ab,a_ind,b_ind = np.intersect1d(b,a,return_indices=True)
# check this also works with datetime

# x = np.array([3,5,7,1,9,8,6,6])
# y = np.array([2,1,5,10,100,6])

# index = np.argsort(x)
# sorted_x = x[index]
# sorted_index = np.searchsorted(sorted_x, y)

# yindex = np.take(index, sorted_index, mode="clip")
# mask = x[yindex] != y
# print result

# asorted = np.argsort(a)
# bpos = np.searchsorted(a[asorted], b)
# indices = asorted[bpos]
