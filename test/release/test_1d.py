import datetime
import hm.api as hm
import xarray
import pytest

def test_set_modeltime():
    hm.set_modeltime(
        datetime.datetime(2003,6,1),
        datetime.datetime(2003,6,30),
        datetime.timedelta(days=1)
    )

def test_set_domain():
    modeltime = hm.set_modeltime(
        datetime.datetime(2003,6,1),
        datetime.datetime(2003,6,30),
        datetime.timedelta(days=1)
    )
    hm.set_domain(
        'test_data/qrparm.soil_HWSD_cont_cosbyNew.nc',
        modeltime,
        'field342',
        is_1d=True,
        xy_dimname='land'
    )
    
def test_open_hmdataarray():
    os.chdir('test/release')
    modeltime = hm.set_modeltime(
        datetime.datetime(2003,6,1),
        datetime.datetime(2003,6,30),
        datetime.timedelta(days=1)
    )
    domain = hm.set_domain(
        'test_data/qrparm.soil_HWSD_cont_cosbyNew.nc',
        modeltime,
        'field342',
        is_1d=True,
        xy_dimname='land'
    )
    # domain.mask
    # domain.starttime
    # domain.endtime
    # domain.n_timestep
    # domain.dt
    da = hm.open_hmdataarray(
        'test_data/SWdown_WFDEI_land_200306.nc',
        'SWdown',
        domain,
        True,
        'land'
    )
    
    da = da.select(
        time=slice(
            datetime.datetime(2003,6,1),
            datetime.datetime(2003,6,10)
        )
    )

# import xarray as xr
# # ds = xr.open_dataset('test_data/AgMERRA_2010_tavg.nc4')
# ds = xr.open_dataset('test_data/SWdown_WFDEI_land_200306.nc', decode_times=False)
# times = ds['time']
# # TODO: cftime.real_datetime to pandas.Timestamp
# # Idea is to manually decode time based on allowable time dimension in dataset
# # this would solve the problem with WFDEI data
# import netCDF4 as nc
# timenum = np.array(
#     nc.num2date(times.values, times.units),
#     dtype='datetime64'
# )
# # convert this to a dataset, then use ds.update()
# dst = xr.Dataset({'tstep' : timenum})
# ds.update(dst)

# import netCDF4 as nc
# ds = nc.Dataset('test_data/SWdown_WFDEI_land_200306.nc')
