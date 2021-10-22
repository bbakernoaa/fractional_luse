from osgeo import gdal
import xarray as xr


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


#def get_source(fname,n_bands=1,in_meta=False):
#    gdal.UseExceptions()
#    ds = gdal.Open(fname)
#    band = ds.GetRasterBand(1)
#    class_ar = band.ReadAsArray()
#    gt = ds.GetGeoTransform()
#    pj = get_projection(fname,in_meta=in_meta)



def digitize_raster(fname,n_bands=17,output_fname='bit_raster.tif',in_meta=False):
    import numpy as np
    from osgeo import gdal
    gdal.UseExceptions()

    # Get data from raster with classifications
    ds = gdal.Open(fname)
    band = ds.GetRasterBand(1)
    class_ar = band.ReadAsArray()
    gt = ds.GetGeoTransform()
    #pj = ds.GetProjection()
    pj = _get_projection(fname,in_meta=in_meta)
    ds = band = None  # close

    # Define the raster values for each class, to relate to each band
    class_ids = (np.arange(n_bands) + 1).tolist()

    # Make a new bit rasters
    drv = gdal.GetDriverByName('GTiff')
    ds = drv.Create(output_fname, class_ar.shape[1], class_ar.shape[0], len(class_ids), gdal.GDT_Byte, ['NBITS=1', 'COMPRESS=LZW', 'INTERLEAVE=BAND'])
    ds.SetGeoTransform(gt)
    ds.SetProjection(pj)
    for bidx in range(ds.RasterCount):
        band = ds.GetRasterBand(bidx + 1)
        # create boolean
        selection = (class_ar == class_ids[bidx])
        band.WriteArray(selection.astype('B'))
    ds = band = None  # save, close

def warp_to_newgrid(src_fname,template,output_fname='nam_viirs_igbp_fluse.tif',percent=False):
    src_ds = gdal.Open(src_fname)
    cpy_ds = gdal.Open(template)
    band = cpy_ds.GetRasterBand(1)
    cpy_mask = (band.ReadAsArray() == band.GetNoDataValue())

    # Result raster, with same resolution and position as the copy raster
    dst_ds = drv.Create(output_fname, cpy_ds.RasterXSize, cpy_ds.RasterYSize,
                        len(class_ids), gdal.GDT_Float32, ['INTERLEAVE=BAND'])
    dst_ds.SetGeoTransform(cpy_ds.GetGeoTransform())
    dst_ds.SetProjection(cpy_ds.GetProjection())

    # Do the same as gdalwarp -r average; this might take a while to finish
    gdal.ReprojectImage(src_ds, dst_ds, None, None, gdal.GRA_Average)
    fac = 1.
    if percent:
        fac = 100.
    # Convert all fractions to percent, and apply the same NODATA mask from the copy raster
    NODATA = 0.
    for bidx in range(dst_ds.RasterCount):
        band = dst_ds.GetRasterBand(bidx + 1)
        ar = band.ReadAsArray() * fac
        ar[cpy_mask] = NODATA
        band.WriteArray(ar)
        band.SetNoDataValue(NODATA)

    # Save and close all rasters
    src_ds = cpy_ds = dst_ds = band = None


def reproject_dataset_digitized(dataset,target_geotransform,xsize,ysize,to_wkt=None,method=gdal.GRA_Bilinear,levels=17):
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
    return data


def reproject_dataset(dataset,target_geotransform,xsize,ysize, to_wkt=None,method=gdal.GRA_Bilinear,in_meta=False):
    """
    A function to reproject and resample a GDAL dataset from within
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
    from osgeo import osr
    source_wkt = _get_projection(dataset,in_meta=in_meta)
    # Up to here, all  the projection have been defined, as well as a
    # transformation from the from to the  to :)
    # We now open the dataset
    g = gdal.Open ( dataset )
    # Get the Geotransform vector
    geo_t = g.GetGeoTransform ()
    # Now, we create an in-memory raster
    mem_drv = gdal.GetDriverByName( 'MEM' )
    # The size of the raster is given the new projection and pixel spacing
    # Using the values we calculated above. Also, setting it to store one band
    # and to use Float32 data type.
    dest = mem_drv.Create('', xsize, ysize, 1, gdal.GDT_Float32)
    # Set the geotransform
    dest.SetGeoTransform( target_geotransform )
    dest.SetProjection ( to_wkt )
    # Perform the projection/resampling
    res = gdal.ReprojectImage( g, dest, source_wkt, to_wkt, method )
    data = dest.ReadAsArray()
    dset= g = None
    
    return data
