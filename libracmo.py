import xarray as xr
from libconst import dpm
import pandas as pd
import numpy as np
"""
Library for reading RACMO2.3p2 data
"""

  
def insert_coordinates(da, model="FGRN11"):
    """
    Insert coordinate arrays for geo-referencing
    
    da : DataArray
    returns : DataArray
    """
    if (model == "ZGRN11"):
        with xr.open_dataset("~/Dropbox/data/maskers/Height_latlon_ZGRN11.nc", decode_cf=False) as ds:
            lon2d = ds.lon
            lat2d = ds.lat
    elif (model == "ANT27"):
        with xr.open_dataset("~/Dropbox/data/masker/Height_latlon_ANT27.nc") as ds:
            lon2d = ds.lon
            lat2d = ds.lat
    elif (model == "XPEN055"):
        with xr.open_dataset("~/Dropbox/data/masker/Height_latlon_XPEN055.nc") as ds:
            lon2d = ds.lon
            lat2d = ds.lat
    elif (model=="ANT11"):
        with xr.open_dataset("~/Dropbox/data/masker/Height_latlon_ANT11.nc") as ds:
            lon2d = ds.lon
            lat2d = ds.lat
    else:
        raise NotImplementedError("Model grid not implemented : "+ model)

    ndim = len(da.dims)
    if (ndim == 2):
        da2  = xr.DataArray(da, dims=('lat','lon'), coords={
                                        'lat2d': (('lat', 'lon'), lat2d), 
                                        'lon2d': (('lat', 'lon'), lon2d)})
    elif (ndim == 3):
        try:
            time_coord_name = 'time'
            time_coord = da[time_coord_name]
        except KeyError:
            time_coord_name = list(da.coords)[2]
            time_coord = da[time_coord_name]
        da2  = xr.DataArray(da, dims=(time_coord_name,'lat','lon'), coords={time_coord_name: (time_coord_name, time_coord),
                                        'lat2d': (('lat', 'lon'), lat2d), 
                                        'lon2d': (('lat', 'lon'), lon2d)})
    else:
        raise NotImplementedError("Not implemented, dims = " + str(ndim))
    
    da2.attrs = da.attrs.copy()
    return da2


def read_racmo23p2_clim(varname, **args):
    """
    wrapper for backwards compatability
    TODO: can be removed at some point
    """
    return read_racmo23p2_GRN_clim(varname, **args)
    
    
def read_racmo23p2_AIS_clim(varname, period='1961-1990', clim_type='ymonmean'):
    pass
    
    
def read_racmo23p2_GRN_clim(varname, period='1961-1990', clim_type='ymonmean'):
    """
    read RACMO 2.3p2 11 km Greenland data 
    Data has been processed with CDO into climatologies (ymonmean, ymonstd)
        
    CURRENTLY IMPLEMENTED
        - period 1961-1990 (precomputed)
    
    OPTIONS
        - clim_type : ymonmean / ymonstd
    """
    
    # Precomputed means
    filename = "/gpfs/fs1/work/lvank/data/racmo/racmo23p2_GRN_monthly/{period}/{varname}_{period}_{clim_type}.nc".format(varname=varname, period=period, clim_type=clim_type)
    with xr.open_dataset(filename, decode_times=False) as ds:
        da = ds[varname].squeeze()
            
    print(varname, da.units, period, clim_type)
    da2 = insert_coordinates(da)
    
    # WORKAROUND : scale racmo mass fields by 12 to obtain annual values from monthly data
###    varnames_scale = ('smb','precip','snowmelt','snowfall','runoff','refreeze','subl')
###    
###    if (varname in varnames_scale):
###        print("INFO: RACMO var "+varname+" is multiplied by 12")
###        da2.values = da2.values * 12 # .values multiplication to preserve attributes
    return da2


def read_racmo23p2_GRN_MM(varname, period='1961-1990'):
    """
    Monthly means
    """
    xr.set_options(enable_cftimeindex=True)
    
    # precip.1961-1990.BN_RACMO2.4_FGRN11.MM.nc
    filename = "/gpfs/fs1/work/lvank/data/racmo/racmo23p2_GRN_monthly/{varname}.{period}.BN_RACMO2.4_FGRN11.MM.nc".format(varname=varname, period=period)
    with xr.open_dataset(filename, decode_times=False) as ds:
        da = ds[varname].squeeze().load()
            
    # Fix time index so we can slice, groupby, etc.
    new_time_index = pd.date_range(start='1961-01-01', end='1990-12-31', freq='MS')
    da['time'] = new_time_index
    
    print(varname, da.units, period)
    da2 = insert_coordinates(da)
    return da2
    
    

def area_racmo23p2_GRN():
    with xr.open_dataset("/glade/work/lvank/racmo/racmo23p2_GRN_monthly/FGRN11_Masks.nc", decode_cf=False) as ds:
        racmo_area = ds['gridarea']
    return racmo_area

# def read_racmo23p2_ann(varname):
#     """
#    1961-1990 mean
#    """
#     da = read_racmo23p2_clim(varname)

def mask_data(var_racmo,model,masker):
    """Masks data
    """
    if model == 'ANT27' and masker != 'none':
        with xr.open_dataset("~/Dropbox/data/masker/Height_latlon_ANT27.nc") as ds:
                mask2d = ds[masker]
                
        var_racmo = np.where(mask2d > 0,var_racmo,np.NaN)
    elif model == 'ANT11' and masker != 'none':
        with xr.open_dataset("~/Dropbox/data/masker/Height_latlon_ANT11.nc") as ds:
            mask2d = ds[masker]
        var_racmo = np.where(mask2d >0,var_racmo,np.NaN)
    elif model == 'XPEN055' and masker != 'none':
        with xr.open_dataset("~/Dropbox/data/masker/Height_latlon_XPEN055.nc") as ds:
            mask2d = ds[masker]
        var_racmo = np.where(mask2d >0.5,var_racmo,np.NaN)   
    
    return var_racmo
