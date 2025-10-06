import os
import warnings

# Use rasterio for reading/writing GeoTIFFs and handling georeferencing info
try:
    import rasterio
    from rasterio.enums import Compression
except ImportError:
    warnings.warn("Rasterio not installed. GeoTIFF functionality is disabled.")
    rasterio = None

def imread_single_band(path):
    """Reads a single band from GeoTIFF using rasterio."""
    with rasterio.open(path) as src:
        # Read the first band and the metadata (profile)
        return src.read(1), src.profile

def compress_geotiff(path_in, path_out):
    """
    // src/+lut/compress_geotiff.m
    Compresses a GeoTIFF using PackBits compression while retaining metadata.
    """
    if rasterio is None:
        warnings.warn("Cannot compress: rasterio library is required for GeoTIFF operations.")
        return

    try:
        im, profile = imread_single_band(path_in)
    except Exception as e:
        warnings.warn(f"Error reading input image {path_in}: {e}. Skipping compression.")
        return

    # MATLAB: ver_out/toolboxes logic is replaced by a try-except block on library
    
    # MATLAB: comp.Compression = Tiff.Compression.PackBits;
    # Rationale: PackBits is lossless and fast. In rasterio, we update the profile.
    profile.update(compress=Compression.packbits)
    
    # The MATLAB GeoTIFFTags/GeoKeyDirectoryTag/SpatialRef logic is implicitly handled 
    # by retaining the rasterio profile.
    
    # MATLAB: geotiffwrite(path_out, im, R, 'GeoKeyDirectoryTag', key, 'TiffTags', comp)
    try:
        with rasterio.open(path_out, 'w', **profile) as dst:
            dst.write(im, 1)
            print(f"Compressed GeoTIFF saved to {path_out} using PackBits.")
            
    except Exception as e:
        warnings.warn(f"Failed to write compressed GeoTIFF: {e}. Output not georeferenced/compressed.")
        print(f"copying georeference tags from input image {path_in}")
        warnings.warn(f"Mapping Toolbox equivalent (rasterio) failed. Output GeoTIFF may be corrupted.")