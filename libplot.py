
import libtimeseries
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
def define_global_map(fig, sps=None):
    """
    creates and returns GeoAxes in Orthographic projection,
    
    """
    import matplotlib
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    
    # Plate Carree
    myproj3 = ccrs.PlateCarree()
    
    #myproj3 = ccrs.Mollweide()


    if (sps == None):
        ax = plt.axes(projection=myproj3)
    elif (type(sps) == matplotlib.gridspec.SubplotSpec): 
        ax = fig.add_subplot(sps, projection=myproj3)
    else: 
        raise TypeError("sps not of type matplotlib.gridspec.SubplotSpec")
    
    #ax.coastlines(resolution='50m') # high res
    ax.coastlines()
        
    ax.set_global()
    #ax.set_extent([-54.5, -27, 59, 84], crs=ccrs.PlateCarree())
    
    gridlines = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, alpha=0.5) 
    
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.0)
    return ax



def get_greenland_shape():
    """
    Get greenland shape, for drawing coastline of Greenland (and no others)
    
    based off
    https://gis.stackexchange.com/questions/88209/python-mapping-in-matplotlib-cartopy-color-one-country
    
    """
    import cartopy.io.shapereader as shpreader

    shpfilename = shpreader.natural_earth(resolution='50m',
                                          category='cultural',
                                          name='admin_0_countries')
    reader = shpreader.Reader(shpfilename)
    countries = reader.records()

    for country in countries:
        if country.attributes['NAME_LONG'] == 'Greenland':
            my_shape = country.geometry
    return my_shape


def define_greenland_map(fig, sps=None, extent=[-54.5, -27, 59, 84], hires_coast=False,
                        central_lon=-40, central_lat=72):
    """
    creates and returns GeoAxes in Orthographic projection,
    constrained to the Greater Greenland region
    
    hires_coast : True, False, or '50m_GREENLAND'
    """
    import matplotlib
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    
    # Orthographic 
    # sometimes doesn't work in Cartopy 0.16 because of this bug https://github.com/SciTools/cartopy/issues/1142
    myproj3 = ccrs.Orthographic(central_longitude=central_lon, central_latitude=central_lat) 
    
    # Stereographic
    # suffers from the problem that Missing Values are shown as 0??? 
    #myproj3 = ccrs.Stereographic(central_longitude=-40, central_latitude=72)

    if (sps == None):
        ax = plt.axes(projection=myproj3)
    elif (type(sps) == matplotlib.gridspec.SubplotSpec): 
        ax = fig.add_subplot(sps, projection=myproj3)
    else: 
        raise TypeError("sps not of type matplotlib.gridspec.SubplotSpec")
    

    if (hires_coast == '50m_GREENLAND'):
        # The following bit only adds a coastline for Greenland
        # without all the other crap (Ellesmere, Iceland)
        gris_shape = get_greenland_shape()
        ax.add_geometries(gris_shape, ccrs.PlateCarree(),
                  edgecolor='black', facecolor='white', zorder=0)
        
    elif (hires_coast == '50m_GREENLAND_FG'):
        # Same as above, but putting the coastline at the foreground
        gris_shape = get_greenland_shape()
        ax.add_geometries(gris_shape, ccrs.PlateCarree(),
                  edgecolor='black', facecolor='none', alpha=1., zorder=2)
        
    elif (hires_coast):
        ax.coastlines(resolution='50m') # high res        
    else:
        ax.coastlines() # low res

    ax.set_extent(extent, crs=ccrs.PlateCarree())
    
    #gridlines = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, alpha=0.5) 
    
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.0)
        
    return ax



def define_greater_greenland_map(fig, sps=None, hires_coast=False):
    """
    creates and returns GeoAxes in Orthographic projection,
    constrained to the Greater Greenland region
    
    """    
    return define_greenland_map(fig, sps=sps, extent=[-60, -10, 55, 78], hires_coast=hires_coast)
      
def define_north_greenland_map(fig, sps=None, hires_coast=False):
    """
    detail map of north Greenland
    initially created to highlight grid imprinting of the EC method in our CESM2 smb evaluation paper
    """    
    return define_greenland_map(fig, sps=sps, extent=[-60, -18, 76, 83.5], hires_coast=hires_coast,
                               central_lon=-50, central_lat=79)
    
    
def define_polar_map(fig, lat_min, lat_max, projection, sps=None):
    """
    Internal helper function for both north/south polar map
    """
    from matplotlib.path import Path
    import numpy as np
    from cartopy.io.shapereader import Reader
    from cartopy.feature import ShapelyFeature
    if (sps == None):
        ax = plt.axes(projection=projection)
        
    elif (type(sps) == matplotlib.gridspec.SubplotSpec): 
        ax = fig.add_subplot(sps, projection=projection)
        
    else: 
        raise TypeError("sps not of type matplotlib.gridspec.SubplotSpec")
    
      
    ax.coastlines(resolution='50m')
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'antarctic_ice_shelves_lines', '50m', edgecolor='k',facecolor='none'))
      
    
    # From example: http://scitools.org.uk/cartopy/docs/latest/examples/always_circular_stereo.html
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = Path(verts * radius + center)

    #ax.set_extent([-180, 180, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.set_boundary(circle, transform=ax.transAxes)
    
    return ax

def define_shapefile(ax):

    
    
    return ax

def define_north_polar_map(fig, lat_min = 40., sps = None):
    """
    creates and returns GeoAxes in NorthPolarStereo projection, constrained to 55-90 N
    
    when plotting on a GridSpec, the SubplotSpec can be passed as an argument
    https://matplotlib.org/api/_as_gen/matplotlib.gridspec.GridSpec.html
    
    e.g.
        gs = gridspec.GridSpec(1, 2, figure=fig)
        ax = define_north_polar_map(fig, gs[0,0])
    
    fig     : Matplotlib figure
    lat_min : minimum latitude still shown
    sps     : subplotspec instance
    """
    return define_polar_map(fig, lat_min=lat_min, lat_max=90., projection=ccrs.NorthPolarStereo(), sps=sps)
    



def define_south_polar_map(fig, lat_max = -55., sps = None):
    """
    creates and returns GeoAxes in SouthPolarStereo projection, constrained to 55-90 S
    
    when plotting on a GridSpec, the SubplotSpec can be passed as an argument
    https://matplotlib.org/api/_as_gen/matplotlib.gridspec.GridSpec.html
    
    e.g.
        gs = gridspec.GridSpec(1, 2, figure=fig)
        ax = define_north_polar_map(fig, gs[0,0])
    
    fig     : Matplotlib figure
    lat_min : minimum latitude still shown
    sps     : subplotspec instance
    """
    return define_polar_map(fig, lat_min=-70, lat_max=lat_max, projection=ccrs.SouthPolarStereo(), sps=sps)


def define_peninsula_map(fig, lat_max = -60., sps = None):
    """
    creates and returns GeoAxes in SouthPolarStereo projection, constrained to 55-90 S
    
    when plotting on a GridSpec, the SubplotSpec can be passed as an argument
    https://matplotlib.org/api/_as_gen/matplotlib.gridspec.GridSpec.html
    
    e.g.
        gs = gridspec.GridSpec(1, 2, figure=fig)
        ax = define_north_polar_map(fig, gs[0,0])
    
    fig     : Matplotlib figure
    lat_min : minimum latitude still shown
    sps     : subplotspec instance
    """

    return define_polar_map(fig, lat_min=-70, lat_max=-60, projection=ccrs.SouthPolarStereo(), sps=sps)


def add_cyclic_point(xarray_obj, dim, period=None):
    """
    to prevent a gap at 0 longitude in Cartopy plot, add one extra longitude 
    to xarray object, the so-called cyclic point
    
    This code was adapted from https://github.com/pydata/xarray/issues/1005
    
    usage: 
        erai_jja = libplot.add_cyclic_point(erai_jja, dim='lon')
    """
    import xarray as xr
    
    if period is None:
        #period = xarray_obj.sizes[dim] / xarray_obj.coords[dim][:2].diff(dim).item()
        period = xarray_obj.sizes[dim] * xarray_obj.coords[dim][:2].diff(dim).item() # LvK
    first_point = xarray_obj.isel({dim: slice(1)})
    first_point.coords[dim] = first_point.coords[dim]+period
    return xr.concat([xarray_obj, first_point], dim=dim)

def get_levels(VARname):
    """
    sets levels, gives units on given RACMO variable name
    """
    masker = 'none'
    if VARname=="precip" or VARname=="smb" or VARname=="snowfall" or VARname=="rainfall":    		
        if VARname=="rainfall":
            levels = [0,1,2,5,10,20,50,100,200,500]
            levelsdiff = [-100,-50,-20,-10,-5,-2,-1,1,2,5,10,20,50,100]
        else:
            levels = [0,2, 5,10,20, 50, 100, 200, 300, 500, 700, 1000, 2000,5000]
            levelsdiff = [-500,-250,-100,-50,-20,-10,-5,-2,-1,1,2,5,10,20,50,100,250,500]
        units = "mm w.e. $y^{-1}$"
        cmap = 'rainbow'
        cmap_diff = 'RdBu_r'
        multi = 12
        
    elif VARname=="snowmelt" or VARname=="refreeze" or VARname=="runoff":
        levels=[0,1,2,5,10,20,50,100,200,500,1000]
        levelsdiff=[-200,-100,-50,-20,-10,-5,-2,-1,0,1,2,5,10,20,50,100,200]
        units = "mm w.e. y$^{-1}$"
        cmap = 'rainbow'
        cmap_diff = 'RdBu_r'
        multi = 12
        
    elif VARname=="subl" or VARname=="suds":
        levels=[-500,-200,-100,-50,-20,-10,-5,-2,2,5,10,20,50,100,200,500]
        levelsdiff=[-250,-200,-150,-100,-50,-20,20,50,100,150,200,250]
        units = "mm w.e. $y^{-1}$"
        cmap  = 'RdBu'
        cmap_diff = 'RdBu'
        if VARname=="subl":
            multi = 12
        else: 
            multi = -12
            
        masker = 'mask2d'
        
    elif VARname=="totpore":
        levels = np.arange(0,50,2)
        levelsdiff = np.arange(-10,10,2)
        units = "m"
        cmap = 'rainbow'
        cmap_diff = 'RdBu_r'
        multi = 1
        
    elif VARname=="tskin" or VARname=="t2m" or VARname=="asn17":
        levels = np.arange(230,280,5)
        levelsdiff =(-10,-8,-6,-4,-2,-1,-0.5,0.5,1,2,4,6,8,10)
        units = "K"
        cmap = 'Blues'
        cmap_diff = 'RdBu_r'
        multi = 1 
        masker = 'mask2d'
    elif VARname=="alb":
        levels = np.arange(0.6,0.9,0.05)
        levelsdiff =(-0.3,-0.2,-0.1,-0.05,-0.02,-0.01,0.01,0.02,0.05,0.1,0.2,0.3)
        units = "-"
        cmap = 'Blues'
        cmap_diff = 'RdBu_r'
        multi = 1 
        masker = 'mask2d'   
    elif VARname=="MOA":
        levels = np.arange(0,2,0.25)
        levelsdiff =(-0.3,-0.2,-0.1,-0.05,-0.02,-0.01,0.01,0.02,0.05,0.1,0.2,0.3)
        units = "-"
        cmap = 'viridis'
        cmap_diff = 'RdBu_r'
        multi = 1 
        masker = 'mask2d' 
    else:
        levels = np.arange(0,1000,100)
        levelsdiff = np.arange(-100,100,10)
        units = ""
        cmap = 'viridis'
        cmap_diff = 'RdBu_r'
        multi = 1
        masker = 'none'
        
        print("For this variable levels should still be defined")
        
    return levels,levelsdiff,units,cmap,cmap_diff,multi,masker


