import os,geopandas,pynhd

def set_streams(**kwargs):
    domain        = kwargs.get('domain',        None)
    fname_streams = kwargs.get('fname_streams', None)
    verbose       = kwargs.get('verbose',       False)
    overwrite     = kwargs.get('overwrite',     False)
    if verbose: print('calling set_streams')
    if domain is None or not isinstance(domain,geopandas.GeoDataFrame):
        raise ValueError(f'set_streams missing required argument domain or is not valid geopandas.GeoDataFrame')
    if fname_streams is None:
        raise ValueError(f'set_streams missing required argument fname_streams')
    if not os.path.isfile(fname_streams) or overwrite:
        if verbose: print(f' using pynhd to download NHDPlusHR flowlines - saving to {fname_streams}')
        nhd = pynhd.NHDPlusHR("flowline").bygeom(geom   =domain.total_bounds,
                                                 geo_crs=domain.crs.to_string())
        nhd = nhd.to_crs(domain.crs)
        nhd = geopandas.clip(nhd, domain.geometry.union_all())
        nhd.to_file(fname_streams, driver="GPKG")
    else:
        if verbose: print(f' using existing NHD line file {fname_streams}')