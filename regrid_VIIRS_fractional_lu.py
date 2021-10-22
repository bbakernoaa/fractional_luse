# coding: utf-8
import gdal_tools as gt
import gdal
import xarray as xr
from pyproj import Proj

def _get_projection(fname,in_meta=False, wkt=True):
    from osgeo import gdal,osr
    srs = osr.SpatialReference()
    f = gdal.Open(fname)
    if in_meta:
        wkt_text = f.GetMetadata()['ssm#esri_pe_string']
    else:
        wkt_text = f.GetProjection()
    srs.ImportFromWkt(wkt_text)
    del f
    if wkt:
        return srs.ExportToWkt()
    else:
        return srs.ExportToProj4()

def reproject_dataset_digitized(dataset,target_geotransform,xsize,ysize,to_wkt=None,method=gdal.GRA_Bilinear,levels=20):
    """
    A function to rasterize, reproject and resample a GDAL dataset from within
    Python. The idea here is to reproject from one system to another, as well as
    changing the geotransform. The procedure is slightly long-winded, but
    goes like this:

    1. Set up the two Spatial Reference systems.
    2. Open the original dataset, and get the geotransform
    3. Calculate bounds of new geotransform by projecting the UL corners
    4. Calculate the number of pixels with the new projection & spacing
    5. Create an in-memory raster dataset
    6. Perform the projection
    """
    from numpy import arange
    source_wkt = _get_projection(dataset)
    dset = gdal.Open(dataset)
    band = dset.GetRasterBand(1)
    class_ar = band.ReadAsArray()
    class_ids = (arange(levels) + 1).tolist()
    mem_drv = gdal.GetDriverByName( 'MEM' )
    ds = mem_drv.Create('', class_ar.shape[1], class_ar.shape[0], len(class_ids), gdal.GDT_Byte, ['INTERLEAVE=BAND'])
    ds.SetGeoTransform(dset.GetGeoTransform())
    ds.SetProjection(dset.GetProjection())
    dset = None
    for bidx in range(ds.RasterCount):
        band = ds.GetRasterBand(bidx + 1)
        selection = (class_ar == class_ids[bidx])
        band.WriteArray(selection.astype('B'))
    dst = mem_drv.Create('', xsize, ysize, levels, gdal.GDT_Float32)
    dst.SetGeoTransform(target_geotransform)
    dst.SetProjection(to_wkt)
    gdal.ReprojectImage(ds,dst,None,to_wkt,method)
    ds = None
    data = dst.ReadAsArray()
    dst = None
    
    return data[:,::-1,::-1]

if __name__ == '__main__':
    VIIRS_FILE = 'VIIRS_GST_2012to2019v3L2_20types_30arcsec.nc' # VIIRS 20 category
    levels = 20 # 20 categories
    target_grid_spec = 'GSD_13km_grid_spec.nc'

    output_name = 'GSD_13km_grid_spec_fractional_landuse.nc'

    # first open the target grid and get target geotransform
    target_projection = '+proj=lcc +lon_0=-97.5 +lat_0=38.5 +lat_1=38.5 +lat_2=38.5 +ellps=GRS80'

    target = xr.open_dataset(target_grid_spec)
        
    p = Proj(target_projection)
    target_wkt = p.to_wkt()

    # assume that the upper left corner is at 0,0 and lon is indexed like lon(y,x)
    x,y = p(target.grid_lont,target.grid_latt) # lat and lon names are grid_lont and grid_latt
    x_upper,y_upper = p(target.lon[0,0],target.lat[0,0])
    pixel_size_x = x[0,1] - x[0,0]
    pixel_size_y = y[1,0] - y[0,0]
    width,height = len(q.grid_xt),len(q.grid_yt)
    gdal_transform = (x_upper,pixel_size_x,0,y_upper,0,pixel_size_y)

    data = reproject_dataset_digitized(VIIRS_FILE,gdal_transform,width,height,to_wkt=target_wkt,method=gdal.GRA_Bilinear,levels=levels)

    q['frac_lu'] = (('lu_cat','grid_yt','grid_xt'), data)
    q.frac_lu.attrs['long_name'] = 'Fractional Land Use'
    q.frac_lu.attrs['_FillValue'] = -1
    q.to_netcdf(output_name)
