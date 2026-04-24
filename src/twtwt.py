import math
import os
import sys
import math
import datetime
import geopandas
import hf_hydrodata
import rasterio 
import numpy
import twtnamelist
import rioxarray

def hf_query_nc(*,
    dt_start: datetime.datetime,
    dt_end: datetime.datetime,
    savedir: str,
    verbose: bool = False,
    huc_id: str = None,
    domain: geopandas.GeoDataFrame = None):

    if verbose: print('calling hf_query')
    if not os.path.isdir(savedir): 
        os.makedirs(savedir,exist_ok=True)
    start_date_str = dt_start.strftime('%Y-%m-%d')
    end_date_str   = (dt_end + datetime.timedelta(days=1)).strftime('%Y-%m-%d')
    options_wtd    = {"dataset"             : "conus1_baseline_mod",
                      "temporal_resolution" : "daily",
                      "start_time"          : start_date_str,
                      "end_time"            : end_date_str}
    if huc_id is not None: 
        options_wtd['huc_id'] = huc_id
    else:                  
        options_wtd['grid_bounds'] = _get_parflow_conus1_bbox(domain)
    hf_hydrodata.get_gridded_files(options_wtd,
                                    variables=['water_table_depth'],
                                    filename_template=os.path.join(savedir,"{dataset}_{variable}_{wy}.nc"),
                                    verbose=verbose)

def hf_query(*,
    dt_start: datetime.datetime,
    dt_end: datetime.datetime,
    huc_id: str = None,
    domain: geopandas.GeoDataFrame = None):

    options_wtd    = {"dataset"             : "conus1_baseline_mod",
                      "variable"            : "water_table_depth",
                      "temporal_resolution" : "daily",
                      "start_time"          : dt_start.strftime('%Y-%m-%d'),
                      "end_time"            : (dt_end + datetime.timedelta(days=1)).strftime('%Y-%m-%d')}
    if huc_id is not None: options_wtd['huc_id']      = huc_id
    else:                  options_wtd['grid_bounds'] = _get_parflow_conus1_bbox(domain)
    hf_data = hf_hydrodata.get_gridded_data(options_wtd)
    if hf_data is None: 
        raise Exception(f'hf_query call to hf_hydrodata.get_gridded_data failed - result is None')
    expected_days = (dt_end - dt_start).days + 1
    if hf_data.shape[0] != expected_days:
        raise Exception(f'hf_hydrodata returned data of unexpected time length or invalid structure')
    return hf_data

def set_wtd_get_flag(*,
    dt_start: datetime.datetime,
    dt_end: datetime.datetime,
    dir_wtd: str,
    overwrite: bool = False):
    
    download_flag = False
    idt = dt_start
    while idt <= dt_end:
        fname = os.path.join(dir_wtd,'wtd_'+idt.strftime('%Y%m%d')+'.tiff')
        if not os.path.isfile(fname) or overwrite:
            download_flag = True
            break
        idt += datetime.timedelta(days=1)
    return download_flag

def break_conus1_tiffs(*,
                       domain: geopandas.GeoDataFrame,
                       dt_start: datetime.datetime,
                       dt_end: datetime.datetime,
                       wtd_in_dir: str,
                       wtd_out_dir: str,
                       verbose: bool = False,  
                       overwrite: bool = False,
                       compress: str = 'zstd',
                       level: int = 9):
    
    if verbose: print('calling break_conus1_tiffs')
    os.makedirs(wtd_out_dir,exist_ok=True)
    conus1_proj, _, _, _ = _get_parflow_conus1_grid_info()
    geom = domain.to_crs(conus1_proj).geometry.union_all()
    idt = dt_start 
    while idt <= dt_end:
        fname_in = os.path.join(wtd_in_dir,f'conus1_baseline_mod_water_table_depth_{idt.strftime("%Y%m%d")}.tiff')
        fname_out = os.path.join(wtd_out_dir,'wtd_'+idt.strftime('%Y%m%d')+'.tiff')
        if not os.path.isfile(fname_out) or overwrite:
            if not os.path.isfile(fname_in):
                raise Exception(f'break_conus1_tiffs could not find conus1 file {fname_in}')
            with rioxarray.open_rasterio(fname_in, masked=True, chunks={'band':1,'x':500,'y':500}) as riox_i:                    
                wtd_dt = riox_i.rio.clip(geometries  = [geom],
                                         crs         = riox_i.rio.crs,
                                         all_touched = True, 
                                         drop        = True, 
                                         from_disk   = True)
                wtd_dt.rio.to_raster(fname_out,compress=compress,zstd_level=level)
        idt += datetime.timedelta(days=1)

def download_hydroframe_data(*,
    domain: geopandas.GeoDataFrame,
    dt_start: datetime.datetime,
    dt_end: datetime.datetime,
    dir_wtd: str,
    verbose: bool = False,
    overwrite: bool = False,
    compress: str = 'zstd',
    level: int = 9):

    if verbose: print('calling download_hydroframe_data')
    if not os.path.isdir(dir_wtd): 
        os.makedirs(dir_wtd,exist_ok=True)
    if verbose: 
        print(f' using hf_hydrodata to download parflow water table depth simulations to {dir_wtd}')
    conus1_proj, _, conus1_transform, conus1_shape = _get_parflow_conus1_grid_info()
    domain = geopandas.GeoDataFrame(domain.drop(columns=['geometry']), 
                                    geometry=domain.to_crs(conus1_proj).buffer(distance=1000),              
                                    crs=conus1_proj)
    hf_data = hf_query(dt_start=dt_start, 
                       dt_end=dt_end, 
                       domain=domain)
    hf_conus1grid_temp = numpy.empty(shape=conus1_shape,
                                     dtype=numpy.float64)
    grid_bounds = _get_parflow_conus1_bbox(domain)
    for i in range(hf_data.shape[0]):
        idt = dt_start + datetime.timedelta(days=i)
        fname = os.path.join(dir_wtd,'wtd_'+idt.strftime('%Y%m%d')+'.tiff')
        if not os.path.isfile(fname) or overwrite:
            hf_conus1grid_temp[grid_bounds[1]:grid_bounds[3],
                               grid_bounds[0]:grid_bounds[2]] = hf_data[i,:,:]
            with rasterio.io.MemoryFile() as memfile:
                hf_conus1data = memfile.open(driver    = "GTiff", 
                                             height    = hf_conus1grid_temp.shape[0], 
                                             width     = hf_conus1grid_temp.shape[1], 
                                             crs       = conus1_proj, 
                                             transform = conus1_transform, 
                                             nodata    = numpy.nan, 
                                             count     = 1, 
                                             dtype     = numpy.float64)
                hf_conus1data.write(hf_conus1grid_temp,1)
                wtd_data, wtd_transform = rasterio.mask.mask(dataset    = hf_conus1data, 
                                                            shapes      = [domain.geometry.union_all()], 
                                                            crop        = True, 
                                                            all_touched = True, 
                                                            filled      = True, 
                                                            pad         = True,
                                                            nodata      = numpy.nan)
                wtd_meta = hf_conus1data.meta
                wtd_meta.update({"driver"   : "GTiff",
                                "height"    : wtd_data.shape[1],
                                "width"     : wtd_data.shape[2],
                                "transform" : wtd_transform, 
                                "nodata"    : numpy.nan,
                                "compress"  : compress,
                                "zlevel"    : level})
                with rasterio.open(fname,'w',**wtd_meta) as wtd_dataset:
                    wtd_dataset.write(wtd_data[0,:,:],1)
    del hf_data, hf_conus1grid_temp

def _get_latlon_parflow_grid(grid_minx,grid_miny,grid_maxx,grid_maxy):
    """Get latlon bbox from ParFlow CONUS1 grid xy bbox"""
    latlon_bounds = hf_hydrodata.to_latlon("conus1", *[grid_minx,grid_miny,grid_maxx,grid_maxy])
    lon_min = latlon_bounds[1]
    lat_min = latlon_bounds[0]
    lon_max = latlon_bounds[3]
    lat_max = latlon_bounds[2]
    return lon_min,lat_min,lon_max,lat_max

def _get_parflow_conus1_bbox(domain:geopandas.GeoDataFrame):
    latlon_tbounds = domain.to_crs(epsg=4326).total_bounds
    conus1grid_minx, conus1grid_miny = hf_hydrodata.from_latlon("conus1", latlon_tbounds[1], latlon_tbounds[0])
    conus1grid_maxx, conus1grid_maxy = hf_hydrodata.from_latlon("conus1", latlon_tbounds[3], latlon_tbounds[2])
    conus1grid_minx, conus1grid_miny = math.floor(conus1grid_minx), math.floor(conus1grid_miny)
    conus1grid_maxx, conus1grid_maxy = math.ceil(conus1grid_maxx),  math.ceil(conus1grid_maxy)
    return tuple([conus1grid_minx, conus1grid_miny, conus1grid_maxx, conus1grid_maxy])

def _set_parflow_conus2_bbox(namelist:twtnamelist.Namelist):
    """Set domain ParFlow bounding box"""
    if namelist.options.verbose: print('calling _get_parflow_conus2_bbox')
    if not os.path.isfile(namelist.fnames.domain): sys.exit('ERROR could not find '+namelist.fnames.domain)
    conus2grid_minx, conus2grid_miny = hf_hydrodata.from_latlon("conus2", namelist.bbox_domain.lat_min, namelist.bbox_domain.lon_min)
    conus2grid_maxx, conus2grid_maxy = hf_hydrodata.from_latlon("conus2", namelist.bbox_domain.lat_max, namelist.bbox_domain.lon_max)
    conus2grid_minx, conus2grid_miny = math.floor(conus2grid_minx), math.floor(conus2grid_miny)
    conus2grid_maxx, conus2grid_maxy = math.ceil(conus2grid_maxx),  math.ceil(conus2grid_maxy)
    return conus2grid_minx,conus2grid_miny,conus2grid_maxx,conus2grid_maxy

def _get_parflow_conus1_grid_info():
    """Parflow CONUS1 grid info - see https://hf-hydrodata.readthedocs.io/en/latest/available_grids.html"""
    conus1_proj      = '+proj=lcc +lat_1=33 +lat_2=45 +lon_0=-96.0 +lat_0=39 +a=6378137.0 +b=6356752.31'
    conus1_spatext   = tuple([-121.47939483437318, 31.651836025255015, -76.09875469594509, 50.49802132270979])
    conus1_transform = rasterio.transform.Affine(1000.0,0.0,-1885055.4995,0.0,1000.0,-604957.0654)
    conus1_shape     = (1888,3342)
    return conus1_proj,conus1_spatext,conus1_transform,conus1_shape

def _get_parflow_conus2_grid_info():
    """Parflow CONUS2 grid info - see https://hf-hydrodata.readthedocs.io/en/latest/available_grids.html"""
    conus2_proj      = '+proj=lcc +lat_1=30 +lat_2=60 +lon_0=-97.0 +lat_0=40.0000076294444 +a=6370000.0 +b=6370000'
    conus2_spatext   = tuple([-126.88755692881833, 21.8170599154073, -64.7677149695924, 53.20274381640737])
    conus2_transform = rasterio.transform.Affine(1000.0,0.0,-2208000.30881173,0.0,1000.0,-1668999.65483222)
    conus2_shape     = (3256,4442)
    return conus2_proj,conus2_spatext,conus2_transform,conus2_shape

def get_conus1_tiffs(*,
    dt_start: datetime.datetime,
    dt_end: datetime.datetime,
    savedir: str,
    verbose: bool = False):
    
    if verbose: print('calling get_conus1_tiffs')
    if not os.path.isdir(savedir): 
        os.makedirs(savedir,exist_ok=True)
    options_wtd = {"dataset"             : "conus1_baseline_mod",
                   "temporal_resolution" : "daily"}
    idatetime = dt_start
    while idatetime <= dt_end:
        if verbose: 
            print(f' working on {idatetime.strftime('%Y-%m-%d')}')
        options_wtd["start_time"] = idatetime.strftime('%Y-%m-%d')
        options_wtd["end_time"]   = idatetime.strftime('%Y-%m-%d')
        fname = os.path.join(savedir,
                             f'conus1_baseline_mod_water_table_depth_{idatetime.strftime('%Y%m%d')}.tiff')
        if not os.path.isfile(fname):
                try:
                    hf_hydrodata.get_gridded_files(options_wtd,
                                                variables=['water_table_depth'],
                                                filename_template=os.path.join(savedir,"{dataset}_{variable}_{ymd}.tiff"),
                                                verbose=True)
                    break
                except Exception as e:
                    sys.stderr.write(f'ERROR: hf_hydrodata.get_gridded_files call failed for date {idatetime.strftime('%Y-%m-%d')} with error: {str(e)}\n')
        idatetime += datetime.timedelta(days=1)