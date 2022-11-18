import xarray as xr
import numpy as np

from libconst import dpm
    
    
def insert_2d_coordinates(ds):
    """
    insert 2d latitude / longitude arrays for pcolormesh plotting
    
    Note: this solution suffers from : 
    https://bairdlangenbrunner.github.io/python-for-climate-scientists/matplotlib/pcolormesh-grid-fix.html
    """
    lon = ds.lon
    lat = ds.lat
    lon2d, lat2d = np.meshgrid(lon,lat)
    #print(lon2d.shape)
    #print(len(lon), len(lat))
    ds.coords['lon2d'] = xr.DataArray(lon2d, dims=('lat', 'lon'))
    ds.coords['lat2d'] = xr.DataArray(lat2d, dims=('lat', 'lon'))   
 

def convert_massflux_units(var):
    
    nt = len(var); assert nt%12 == 0. # monthly data
    
#     print(var.units)
    
    if (var.units == 'm/s'):
        # convert m/s to total mm per month
        print("INFO: converting units from m/s to mm/s for variable "+var.name)
        var.values *= 1000.
        var.attrs['units'] = 'mm/s'
        
    if (var.units == 'mm/s' or var.units == 'kg/m2/s'):
        # convert mm/s to total mm per month
        print("INFO: converting units from mm/s to mm for variable "+var.name)
        ndim = len(var.shape)
        if (ndim == 3):
            var.values *= np.tile(dpm, nt//12)[:,np.newaxis,np.newaxis] * 86400
        elif (ndim == 2):
            var.values *= np.tile(dpm, nt//12)[:,np.newaxis] * 86400
        else:
            raise NotImplementedError("ERROR: var dimension = "+str(ndim))
            
        var.attrs['units'] = 'mm'

        
def monthly_to_clim(field,  clim_type='mean', groupby='time.month'):
    """
    given a variable @ monthly frequency, produce climatology
    """
    if (clim_type == 'mean'):
        out = field.groupby(groupby).mean('time')
    elif (clim_type == 'std'):
        out = field.groupby(groupby).std('time')
    elif (clim_type == 'var'):
        out = field.groupby(groupby).var('time')
    else:
        raise(ValueError("unknown clim_type : "+ str(clim_type)))    
    
    # Keep attributes
    out.attrs = field.attrs
    return out


def monthly_to_yearly(var, clim_type='mean'):
    """
    given a variable @ monthly frequency, produce annual values
    """
    
    # Arithmetic mean
    #if (clim_type == 'mean'):
    #    out = var.groupby('time.year').mean('time') 
    #else:
    #    raise(ValueError("unknown clim_type : "+ str(clim_type)))    
    mymask = var.to_masked_array()
        
        
    ## MEAN WEIGHTED BY NUMBER OF DAYS (ASSUMING 365 DAY CALENDAR)
    def mul_with_dpm(x):
        return (x * dpm[:,np.newaxis,np.newaxis]).sum(dim='time') / dpm.sum()

    # Keep attributes
    #out.attrs = var.attrs
    
    if (clim_type == 'mean'):
        out = var.groupby('time.year').apply(mul_with_dpm)
        #out.where(mymask)
        out.attrs = var.attrs
        return out
    else:
        raise(ValueError("unknown clim_type : "+ str(clim_type)))    
        
    
   
    
def read_monthly_from_timeseries(filename, varname, ys=1980, ye=1999, convert_units=True, insert_2d_latlon = True, verbose=True):
    """
    read data from a single CESM monthly timeseries file (CLM / CAM)
    
    filename  : string
    varname   : string
    ys        : starting year
    ye        : ending year
    """
    nyear = ye-ys+1
    
    if (type(filename) == list or type(filename) == tuple):
        #ds = xr.open_mfdataset(filename)
        ds = xr.open_mfdataset(filename,combine='by_coords')
    else:
        ds = xr.open_dataset(filename)
        
    #print (filename)
        
    # shifting time doesnt seem to work for deviating time units
    #assert(type(ds.time.values[0]) == np.datetime64) ## STOPPED WORKING SINCE XARRAY 0.11.1
    
    # shift time coordinate by one day
    # this is needed because CESM puts the time stamp at the very end of the
    # averaging period, see https://bb.cgd.ucar.edu/external-tools-dont-cesm-time-coordinate
    from datetime import timedelta
    time_index_shifted  = ds.time.get_index('time') - timedelta(days=1)
    ds['time'] = time_index_shifted
    #print(time_index_shifted)
    
    #print(ds.time)
    my_units = ds[varname].attrs['units']
    
    if (verbose):
        print(varname, my_units, 'ys, ye, nyear', ys, ye, nyear)

    if (insert_2d_latlon): 
        insert_2d_coordinates(ds)
    
    
    ### SELECT VARIABLE OF INTEREST
    var = ds[varname].sel(time=slice('{:04d}'.format(ys), '{:04d}'.format(ye))).squeeze().load()
    assert len(var) == nyear * 12, 'len(var) = %d' % len(var)
    
    if (convert_units):
        convert_massflux_units(var)
    else:
        var.attrs['units'] = my_units
        
    ds.close()
    return var


def read_clim_from_timeseries(filename, varname, clim_type='mean', ys=1980, ye=1999, convert_units=True, insert_2d_latlon = True):
    """
    read single CESM timeseries file and compute the climatology (12 months)
    
    filename         : string or list of strings
    varname          : string
    clim_type        : type of climatology (mean / std)
    ys               : starting year
    ye               : ending year
    convert_units    : convert mass fluxes to total mass
    insert_2d_latlon : insert 2d coordinate arrays for plotting purposes
    """
    nyear = ye-ys+1
    
    var = read_monthly_from_timeseries(filename, varname=varname, ys=ys, ye=ye, 
                                       convert_units=convert_units, insert_2d_latlon=insert_2d_latlon)
    assert len(var) == nyear * 12 # e.g. 20 years

    return monthly_to_clim(var.squeeze(), clim_type)


def read_yearmean_from_timeseries(filename, varname, ys=1980, ye=1999, convert_units=True, verbose=True, insert_2d_latlon = True):
    """
    read CESM timeseries file and compute yearly means
    
    filename  : string
    varname   : string
    ys        : starting year
    ye        : ending year
    """
    nyear = ye-ys+1
    
    var = read_monthly_from_timeseries(filename, varname, ys, ye, convert_units=convert_units, 
                                       insert_2d_latlon=insert_2d_latlon, verbose=verbose)
    assert len(var) == nyear * 12 # e.g. 20 years
        
    out = monthly_to_yearly(var)
    return out.squeeze()


def read_lon_lat_from_timeseries(filename):
    """
    read CESM longitude latitude info
    """
    if (type(filename) == list or type(filename) == tuple):
        with xr.open_mfdataset(filename) as ds:
            return ds.lon, ds.lat
    else:
        with xr.open_dataset(filename) as ds:
            return ds.lon, ds.lat
        

def month_to_season(var, season, agg_type='mean', arithmetic_mean=False):
    """
    convert monthly variable (leading dimension 12) to seasonal mean
    
    var      :  xarray DataArray
    season   :  string
    agg_type :  aggregation type (mean or sum)
    
    returns: xarray DataArray
    
    """
    assert len(var)==12, 'data is not seasonal? Length = {}'.format(len(var))
    var = var.load()
    
    assert season in ('ANN', 'MAM', 'JJA', 'SON', 'DJF')
    
    if (agg_type == 'mean' and arithmetic_mean):
        ## ARITHMETIC MEAN
        if (season == "ANN"):
            return var.mean(axis=0)
        elif (season == "MAM"):
            return var[[2,3,4]].mean(axis=0)
        elif (season == "JJA"):
            return var[[5,6,7]].mean(axis=0)
        elif (season == "SON"):
            return var[[8,9,10]].mean(axis=0)
        elif (season == "DJF"):
            return var[[0,1,11]].mean(axis=0)
        
    if (agg_type == 'mean' and not arithmetic_mean):
        ## MEAN WEIGHTED BY NUMBER OF DAYS (ASSUMING 365 DAY CALENDAR)
        assert dpm.sum() == 365
        
        wvar = var * dpm[:,np.newaxis,np.newaxis] 
        
        if (season == "ANN"):
            return wvar.sum(axis=0, skipna=False) / 365.
        elif (season == "MAM"):
            return wvar[[2,3,4]].sum(axis=0, skipna=False) / dpm[[2,3,4]].sum()
        elif (season == "JJA"):
            return wvar[[5,6,7]].sum(axis=0, skipna=False)  / dpm[[5,6,7]].sum()
        elif (season == "SON"):
            return wvar[[8,9,10]].sum(axis=0, skipna=False) / dpm[[8,9,10]].sum()
        elif (season == "DJF"):
            return wvar[[0,1,11]].sum(axis=0, skipna=False) / dpm[[0,1,11]].sum()
        
    elif (agg_type == 'sum'):
        ### ONLY TO BE USED FOR SUMMING MASS FLUXES
        ### NO WEIGHTING APPLIED SINCE THESE FIELDS ARE PROBABLY
        ### ALREADY IN TOTAL MM PER MONTH,  WHICH WE SHOULD ASSERT
        assert var.attrs['units'] == 'mm' or var.attrs['units'] == 'mmWE'
        
        if (season == "ANN"):
            return var.sum(axis=0, skipna=False)
        elif (season == "MAM"):
            return var[[2,3,4]].sum(axis=0, skipna=False)
        elif (season == "JJA"):
            return var[[5,6,7]].sum(axis=0, skipna=False)
        elif (season == "SON"):
            return var[[8,9,10]].sum(axis=0, skipna=False)
        elif (season == "DJF"):
            return var[[0,1,11]].sum(axis=0, skipna=False)
        
    else:
        raise ValueError('unknown agg_type: '+agg_type)
    

    
    