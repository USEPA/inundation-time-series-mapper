import os,pygeohydro,geopandas,shapely
    
def set_domain(*,
    fname_domain: str,
    verbose: bool = False,
    overwrite: bool = False,
    domain_hucid: str = None,
    domain_bbox: tuple = None,
    domain_latlon: list = None,
    conus1_domain: str = None):
    if verbose: print('calling set_domain')
    if fname_domain is None: 
        raise ValueError(f'_set_domain missing required argument fname_domain')
    if not os.path.isfile(fname_domain) or overwrite:
        if domain_hucid is not None:
            domain = _set_domain_byhucid(fname_domain=fname_domain, 
                                         domain_hucid=domain_hucid, 
                                         verbose=verbose)
        elif domain_bbox is not None:
            domain = _set_domain_bybbox(fname_domain=fname_domain, 
                                        domain_bbox=domain_bbox, 
                                        verbose=verbose)
        elif domain_latlon is not None:
            domain = _set_domain_bylatlonandhuclvl(fname_domain=fname_domain, 
                                                   domain_latlon=domain_latlon, 
                                                   verbose=verbose)
        else:
            raise Exception(f'_set_domain could not set domain from arguments')
        if conus1_domain is not None:
            conus1_domain = geopandas.read_file(conus1_domain)
            domain = geopandas.clip(gdf=domain,mask=conus1_domain.to_crs(domain.crs))
            domain.to_file(fname_domain, driver='GPKG')
    else:
        print(f' using existing domain {fname_domain}')
        domain = geopandas.read_file(fname_domain)
    return domain
        
def set_domain_buf(*,
    domain: geopandas.GeoDataFrame,
    buf_dist_m: int = 1000,
    fname_domain_buf: str = None,
    verbose: bool = False,
    overwrite: bool = False):
    if verbose: print('calling set_domain_buf')
    if not os.path.isfile(fname_domain_buf) or overwrite:
        if verbose: print(f' creating domain buffer {fname_domain_buf} with buffer distance {buf_dist_m} m')
        domain_buf = geopandas.GeoDataFrame(domain.drop(columns=['geometry']), 
                                                        geometry=domain.to_crs(crs="EPSG:5070").buffer(distance=buf_dist_m),              
                                                        crs="EPSG:5070").to_crs(crs=domain.crs)
        os.makedirs(name=os.path.dirname(fname_domain_buf),exist_ok=True)
        domain_buf.to_file(fname_domain_buf, driver='GPKG')
    else:
        if verbose: print(f' found existing domain buffer file {fname_domain_buf}')
        domain_buf = geopandas.read_file(fname_domain_buf)
    return domain_buf

def _set_domain_byhucid(*,
    fname_domain: str,
    domain_hucid: str,
    verbose: bool = False):
    if verbose: print('_set_domain_byhucid')
    if len(domain_hucid) not in (2,4,6,8,10,12):
        raise ValueError(f'_set_domain_byhucid domain_hucid {domain_hucid} is invalid, must be of len 2, 4, 6, 8, 10, 12')
    colnam = f'huc{len(domain_hucid)}'
    hucs   = pygeohydro.WBD(colnam) 
    domain = hucs.byids(colnam, domain_hucid, return_geom=True)
    domain = domain.drop(columns=[col for col in domain.columns if col not in [colnam,'geometry']]) 
    domain = domain.rename(columns = {colnam : 'domain_id'})
    domain.geometry = domain.geometry.force_2d()
    os.makedirs(name=os.path.dirname(fname_domain),exist_ok=True)
    domain.to_file(fname_domain, driver='GPKG')
    return domain

def _set_domain_bybbox(*,
    fname_domain: str,
    domain_bbox: tuple,
    verbose: bool = False):
    if verbose: print('_set_domain_bybbox')
    geom = shapely.geometry.box(domain_bbox[0],
                                domain_bbox[1],
                                domain_bbox[2],
                                domain_bbox[3])
    domain = geopandas.GeoDataFrame(geometry=[geom], crs="EPSG:4326") # lat/lon
    os.makedirs(name=os.path.dirname(fname_domain),exist_ok=True)
    domain.to_file(fname_domain, driver='GPKG')
    return domain

def _set_domain_bylatlonandhuclvl(*,
    fname_domain: str,
    domain_latlon: list,
    huc_lvl: int = 12,
    verbose: bool = False):
    if verbose: print('_set_domain_bybbox')
    if int(huc_lvl) not in (2,4,6,8,10,12):
        raise ValueError(f'_set_domain_bylatlonandhuclvl huc_lvl {huc_lvl} is invalid, must be 2, 4, 6, 8, 10, 12')
    if not isinstance(domain_latlon,list) or not all(isinstance(v, float) for v in domain_latlon) or len(domain_latlon) != 2:
        raise ValueError(f'_set_domain_bylatlonandhuclvl invalid domain_latlon')
    geom = shapely.geometry.Point(domain_latlon[1],domain_latlon[0]).buffer(0.01)
    colnam = f'huc{huc_lvl}'
    hucs   = pygeohydro.WBD(colnam)
    domain = hucs.bygeom(geom)
    domain = domain.drop(columns=[col for col in domain.columns if col not in [colnam,'geometry']]) 
    domain.geometry = domain.geometry.force_2d()
    os.makedirs(name=os.path.dirname(fname_domain),exist_ok=True)
    domain.to_file(fname_domain, driver='GPKG')
    return domain

def get_conus1_hucs(*, 
    fname_domain: str, 
    fname_domain_hucs: str, 
    huc_lvl: int = 8, 
    verbose: bool = False):
    if verbose: print('calling get_conus1_hucs')
    if not os.path.isfile(fname_domain_hucs):
        fname_wb_full_temp = str(os.path.join(os.path.dirname(fname_domain_hucs),
                                              f'wb_full_huc{str(huc_lvl)}.gpkg'))
        if os.path.isfile(fname_wb_full_temp): 
            wb_full = geopandas.read_file(fname_wb_full_temp)
        else:
            wb_full = pygeohydro.watershed.huc_wb_full(huc_lvl)
            wb_full.to_file(fname_wb_full_temp)
        domain      = geopandas.read_file(fname_domain)
        domain_hucs = geopandas.clip(gdf=wb_full,mask=domain.to_crs(wb_full.crs))
        domain_hucs = domain_hucs[~domain_hucs['states'].str.contains('CN')]
        domain_hucs = domain_hucs[domain_hucs['name'] != 'Lake Michigan']
        domain_hucs = domain_hucs.drop(columns=[col for col in domain_hucs.columns if col not in [f'huc{huc_lvl}','geometry']]) 
        domain_hucs.geometry = domain_hucs.geometry.force_2d()
        domain_hucs.to_file(fname_domain_hucs,driver='GPKG')
    else:
        domain_hucs = geopandas.read_file(fname_wb_full_temp)
    return domain_hucs