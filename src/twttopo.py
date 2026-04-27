import os
import rasterio
import py3dep
import geopandas
import rioxarray
import whitebox
import pynhd
import glob
from osgeo import gdal
from geocube.api.core import make_geocube
import math
import numpy as np
from affine import Affine
from rasterio.warp import transform_bounds
from asyncio import sleep
from py3dep.exceptions import ServiceUnavailableError

async def download_dem(*,
          domain: geopandas.GeoDataFrame,
          fname_dem: str,
          dem_rez: int = None,
          overwrite: bool = False,
          verbose: bool = False,
          compress: str = 'ztsd',
          zlevel:int    = 19,
          n_retries: int = 10,  # number of retry attempts for DEM download in case of failure
          base_wait: int = 60): # seconds to wait before retrying after a failure, will be multiplied by 2^attempt for exponential backoff
    
    if verbose: print(f'calling download_dem')
    if not os.path.isfile(fname_dem) or overwrite:
        avail = py3dep.check_3dep_availability(bbox=tuple(domain.total_bounds),
                                               crs=domain.crs)
        vals  = list()
        for k in avail.keys():
            try: vals.append(int(k.replace('m','')))
            except: pass
        if len(vals) == 0:
            raise ValueError(f'download_dem 3dep could not find dem resolution')
        rez = min(vals)
        if dem_rez in vals: rez = dem_rez # override with user input if available
        if verbose:
            print(f' downloading {str(rez)}m DEM to {fname_dem}')
        for attempt in range(1, n_retries + 1):
            wait = base_wait * (2 ** attempt)
            try:
                dem = py3dep.get_dem(geometry   = domain.geometry.union_all(),
                                     resolution = rez,
                                     crs        = domain.crs)
                dem.rio.to_raster(fname_dem, 
                                driver="GTiff", 
                                compress=compress, 
                                zlevel=zlevel)
                return
            except ServiceUnavailableError as ex:
                await sleep(wait)
            except TimeoutError as ex:
                if verbose:
                    print(f" DEM download timed out (attempt {attempt}/{n_retries}) at {rez} m")
                if attempt < n_retries:
                    await sleep(wait) 
                else:
                    raise Exception(f"DEM download failed after {n_retries} attempts at {rez} m.") from ex
    else:
        if verbose: print(f' found existing dem {fname_dem}')

def break_dem(*,
    fname_dem_parent: str,
    fname_dem_child: str,
    fname_boundary: str,
    overwrite: bool = False,
    verbose: bool = False):

    if verbose: print('calling break_dem')
    if not os.path.isfile(fname_dem_child) or overwrite:
        if verbose: print(f' creating {fname_dem_child} from {fname_dem_parent}')
        gdal.UseExceptions()
        src_ds = gdal.Open(fname_dem_parent, gdal.GA_ReadOnly)
        if src_ds is None: raise RuntimeError(f"Failed to open VRT: {fname_dem_parent}")
        warp_kwargs = {
            "format": "GTiff",
            "cutlineDSName": fname_boundary,      
            "cropToCutline": True,                
            "resampleAlg": "near",             
            "multithread": True,
            "warpOptions": ["NUM_THREADS=ALL_CPUS"],
            "warpMemoryLimit": 512 * 1024 * 1024,
            "creationOptions": [
                "TILED=YES",
                "COMPRESS=LZW",
                "BIGTIFF=IF_SAFER",
            ],
        }
        try:
            nodata = src_ds.GetRasterBand(1).GetNoDataValue()
            warp_kwargs["dstNodata"] = nodata
        except Exception:
            pass
        out_ds = gdal.Warp(fname_dem_child, src_ds, **warp_kwargs)
        if out_ds is None: raise RuntimeError("GDAL Warp failed")
        out_ds = None
        src_ds = None
    else:
        if verbose: print(f' found {fname_dem_child}')

def breach_dem(*,
    fname_dem_breached: str,
    fname_dem: str,
    verbose: bool = False,
    overwrite: bool = False):

    if verbose: print('calling breach_dem')
    if not os.path.isfile(fname_dem_breached) or overwrite:
        if verbose: 
            print(f' using whitebox to breach dem and writing to {fname_dem_breached}')
        wbt = whitebox.WhiteboxTools()
        wbt.set_compress_rasters(True) 
        fname_filled = fname_dem.replace('.tif','_fill.tif')
        wbt.breach_single_cell_pits(
            dem=fname_dem, 
            output=fname_filled, 
        )
        wbt.breach_depressions_least_cost(
            dem=fname_filled,
            output=fname_dem_breached,
            dist=int(1000),
        )
        os.remove(fname_filled)
    else:
        if verbose: print(f' found existing breached dem {fname_dem_breached}')

def set_flow_acc(*,
    fname_dem_breached: str,
    fname_facc_ncells: str,
    fname_facc_sca: str,
    verbose: bool = False,
    overwrite: bool = False):

    if verbose: print(f'calling set_flow_acc')
    if not os.path.isfile(fname_facc_ncells) or overwrite:
        if verbose: print(f' using whitebox to calculate flow accumulation (n cells), writing to {fname_facc_ncells}')
        wbt = whitebox.WhiteboxTools()
        wbt.set_compress_rasters(True) 
        wbt.d_inf_flow_accumulation(i        = fname_dem_breached,
                                    output   = fname_facc_ncells,
                                    out_type = 'cells',
                                    log      = False)
    else:
        if verbose: print(f' using existing flow accumulation (ncells) file {fname_facc_ncells}')
    if not os.path.isfile(fname_facc_sca) or overwrite:
        if verbose: print(f' using whitebox to calculate flow accumulation (sca), writing to {fname_facc_sca}')
        wbt = whitebox.WhiteboxTools()
        wbt.set_compress_rasters(True) 
        wbt.d_inf_flow_accumulation(i        = fname_dem_breached,
                                    output   = fname_facc_sca,
                                    out_type = 'sca',
                                    log      = False)

    else:
        if verbose: print(f' using existing flow accumulation (sca) file {fname_facc_sca}')

def calc_stream_mask(*,
    verbose: bool = False,
    overwrite: bool = False,
    fname_facc_ncells: str = None,
    fname_facc_sca: str = None,
    facc_threshold_ncells: float = None,
    facc_threshold_sca: float = None,
    fname_strm_mask: str = None):

    if verbose: print('calling calc_stream_mask')
    if not os.path.isfile(fname_strm_mask) or overwrite:
        if os.path.isfile(fname_facc_ncells):
            if verbose: print(f' setting stream mask using fname_facc_ncells {fname_facc_ncells} and facc_threshold_ncells {facc_threshold_ncells}')
            wbt = whitebox.WhiteboxTools()
            wbt.set_compress_rasters(True) 
            wbt.extract_streams(flow_accum      = fname_facc_ncells,
                                output          = fname_strm_mask,
                                threshold       = facc_threshold_ncells,
                                zero_background = True)
        elif os.path.isfile(fname_facc_sca):
            if verbose: print(f' setting stream mask using fname_facc_sca {fname_facc_sca} and facc_threshold_sca {facc_threshold_sca}')
            wbt = whitebox.WhiteboxTools()
            wbt.set_compress_rasters(True) 
            wbt.extract_streams(flow_accum      = fname_facc_sca,
                                output          = fname_strm_mask,
                                threshold       = facc_threshold_sca,
                                zero_background = True)
        else:
            raise Exception(f'calc_stream_mask did not find valid flow accumulation file fname_facc_ncells {fname_facc_ncells} or fname_facc_sca {fname_facc_sca}')
    else:
        if verbose: print(f' using existing stream mask {fname_strm_mask}')

def calc_slope(*,
    fname_dem_breached: str,
    fname_slope: str,
    verbose: bool = False,
    overwrite: bool = False):

    if verbose: print(f'calling calc_slope')
    if not os.path.isfile(fname_slope) or overwrite:
        wbt = whitebox.WhiteboxTools()
        wbt.set_compress_rasters(True) 
        wbt.slope(dem    = fname_dem_breached,
                  output = fname_slope,
                  units  = 'degrees')
    else:
        if verbose: print(f' using existing slope file {fname_slope}')

def calc_twi(*,
    fname_facc_sca: str,
    fname_slope: str,
    fname_twi: str,
    verbose: bool = False,
    overwrite: bool = False):

    if verbose: print('calling calc_twi')
    if not os.path.isfile(fname_twi) or overwrite:
        if verbose: print(f' using whitebox workflows to calculate {fname_twi}')
        wbt = whitebox.WhiteboxTools()
        wbt.set_compress_rasters(True) 
        wbt.wetness_index(sca    = fname_facc_sca,
                          slope  = fname_slope,
                          output = fname_twi)
    else:
        if verbose: print(f' found existing TWI file {fname_twi}')

def calc_twi_mean(*,
    fname_twi: str,
    fname_twi_mean: str,
    wtd_raw_dir: str,
    verbose: bool = False,
    overwrite: bool = False,
    compress: str = 'ztsd',
    zlevel: int = 19):
    """
    Compute the mean TWI per WTD grid cell over the full TWI extent,
    then reproject the result back to the original TWI grid so all TWI cells
    are present in the output (no cut-off).

    Required kwargs:
      - fname_twi: path to input TWI raster
      - fname_twi_mean: path to output raster
      - wtd_raw_dir: directory containing one or more WTD rasters (.tif/.tiff)

    Optional kwargs:
      - verbose: bool
      - overwrite: bool
      - compress: str (e.g., 'deflate' or 'lzw')
    """
    if verbose:
        print('calling calc_twi_mean')

    if not fname_twi or not fname_twi_mean or not wtd_raw_dir:
        raise ValueError('fname_twi, fname_twi_mean, and wtd_raw_dir are required')

    if os.path.isfile(fname_twi_mean) and not overwrite:
        if verbose:
            print(f' found existing mean TWI file {fname_twi_mean}')
        return

    # Pick a representative WTD raster to define CRS/resolution/grid alignment
    wtd_candidates = sorted(
        glob.glob(os.path.join(wtd_raw_dir, '*.tif')) +
        glob.glob(os.path.join(wtd_raw_dir, '*.tiff'))
    )
    if not wtd_candidates:
        raise ValueError(f'calc_twi_mean could not locate any .tif/.tiff files in {wtd_raw_dir}')
    fname_example_wtd_raw = wtd_candidates[0]

    if verbose:
        print(f' calculating mean TWI per WTD cell and saving on TWI grid -> {fname_twi_mean}')

    with rioxarray.open_rasterio(fname_twi, masked=True) as riox_ds_twi, \
         rioxarray.open_rasterio(fname_example_wtd_raw, masked=True) as riox_ds_wtd:

        # Validate CRS presence
        if riox_ds_twi.rio.crs is None:
            raise ValueError('Input TWI raster has no CRS defined.')
        if riox_ds_wtd.rio.crs is None:
            raise ValueError('Example WTD raster has no CRS defined.')

        # Build a WTD-grid template that fully covers the TWI extent (prevents cut-off)
        wtd_crs = riox_ds_wtd.rio.crs
        wtd_transform = riox_ds_wtd.rio.transform()

        left_twi, bottom_twi, right_twi, top_twi = riox_ds_twi.rio.bounds()
        twi_crs = riox_ds_twi.rio.crs

        left_wtd, bottom_wtd, right_wtd, top_wtd = transform_bounds(
            twi_crs, wtd_crs, left_twi, bottom_twi, right_twi, top_twi, densify_pts=21
        )

        # Convert those bounds to WTD pixel indices and snap to grid
        inv = ~wtd_transform
        c0, r0 = inv * (left_wtd,  top_wtd)
        c1, r1 = inv * (right_wtd, bottom_wtd)

        col_min = math.floor(min(c0, c1))
        row_min = math.floor(min(r0, r1))
        col_max = math.ceil(max(c0, c1))
        row_max = math.ceil(max(r0, r1))

        width = int(col_max - col_min)
        height = int(row_max - row_min)

        if width <= 0 or height <= 0:
            raise ValueError('Computed WTD template has non-positive dimensions. Check input extents/CRS.')

        template_transform = wtd_transform * Affine.translation(col_min, row_min)

        # Aggregate TWI -> WTD grid using average resampling on the full-coverage template.
        # IMPORTANT: Do not pass src_nodata/dst_nodata; rioxarray manages those internally.
        twi_on_wtd = riox_ds_twi.rio.reproject(
            dst_crs=wtd_crs,
            transform=template_transform,
            shape=(height, width),
            resampling=rasterio.enums.Resampling.average
        )

        # Always project back to the TWI grid so the output has exactly the TWI footprint
        twi_mean = twi_on_wtd.rio.reproject_match(
            riox_ds_twi,
            resampling=rasterio.enums.Resampling.nearest
        )

        # Ensure floating dtype for averaged values
        if np.issubdtype(twi_mean.dtype, np.integer):
            twi_mean = twi_mean.astype('float32')
        else:
            # Keep float32 to reduce file size if float64
            try:
                twi_mean = twi_mean.astype('float32')
            except Exception:
                pass

        # Write output
        os.makedirs(os.path.dirname(fname_twi_mean) or '.', exist_ok=True)
        twi_mean.rio.to_raster(
            fname_twi_mean,
            compress=compress,
            zlevel=zlevel,
            tiled=True
        )

        if verbose:
            print(f' wrote {fname_twi_mean}')

def set_domain_mask(*,
    domain: geopandas.GeoDataFrame,
    fname_domain_mask: str,
    fname_dem: str,
    verbose: bool = False,
    overwrite: bool = False,
    compress: str = 'ztsd',
    zlevel: int = 19):

    if verbose: print(f'calling set_domain_mask')
    if not os.path.isfile(fname_domain_mask) or overwrite:
        if verbose: print(f' creating {fname_domain_mask}')
        with rioxarray.open_rasterio(fname_dem) as riox_ds_dem:
            domain = domain.to_crs(riox_ds_dem.rio.crs)
            domain = domain.dissolve()
            domain = domain.drop(columns=[col for col in domain.columns if col not in ['geometry']])
            domain['mask'] = 1
            domain_mask = make_geocube(vector_data=domain,like=riox_ds_dem,measurements=['mask'])
            domain_mask.rio.to_raster(fname_domain_mask,
                                      driver="GTiff",
                                      compress=compress,
                                      zlevel=zlevel,
                                      dtype='uint8',
                                      nodata=0)
    else:
        if verbose: print(f' using existing domain mask {fname_domain_mask}')

def set_streams(*,
    domain: geopandas.GeoDataFrame,
    fname_streams: str,
    verbose: bool = False,
    overwrite: bool = False):
    
    if verbose: print('calling set_streams')
    if not os.path.isfile(fname_streams) or overwrite:
        if verbose: print(f' using pynhd to download NHDPlusHR flowlines - saving to {fname_streams}')
        nhd = pynhd.NHDPlusHR("flowline").bygeom(geom   =domain.total_bounds,
                                                 geo_crs=domain.crs.to_string())
        nhd = nhd.to_crs(domain.crs)
        nhd = geopandas.clip(nhd, domain.geometry.union_all())
        nhd.to_file(fname_streams, driver="GPKG")
    else:
        if verbose: print(f' using existing NHD line file {fname_streams}')






