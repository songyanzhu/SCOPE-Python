import numpy as np
import os
import warnings

# Use rasterio for GeoTIFF writing
try:
    import rasterio
    from rasterio.enums import Compression
except ImportError:
    warnings.warn("Rasterio not installed. GeoTIFF functionality is disabled.")
    rasterio = None

def write_tiff(im, out_path):
    """
    // src/+lut/write_tiff.m
    Writes a numpy array to a GeoTIFF. It mimics MATLAB's behavior of 
    copying an existing GeoTIFF's metadata (Cab.tif) and then overwriting the data.
    """
    if rasterio is None:
        warnings.warn("Cannot write GeoTIFF: rasterio library is required.")
        return

    # MATLAB: input_image_path = '../exercise/images/Cab.tif'; copyfile(input_image_path, out_path)
    input_image_path = os.path.join('..', 'exercise', 'images', 'Cab.tif')
    
    if not os.path.exists(input_image_path):
        warnings.warn(f"Reference image {input_image_path} not found. Cannot copy metadata. Writing basic GeoTIFF.")
        # Create a basic profile if reference is missing
        profile = {'driver': 'GTiff', 'dtype': im.dtype, 'nodata': np.nan, 
                   'width': im.shape[1], 'height': im.shape[0], 'count': 1,
                   'compress': Compression.packbits}
    else:
        # 1. Read metadata from reference
        with rasterio.open(input_image_path) as src:
            profile = src.profile
        
        # 2. Update profile for new image data and compression
        # MATLAB: if getTag(t, 'BitsPerSample') == 32: im = single(im);
        profile.update(dtype=np.float32) # Ensure 32-bit float
        profile.update(compress=Compression.packbits)
        profile.update(nodata=np.nan)
        profile.update(width=im.shape[1], height=im.shape[0])
        
    # 3. Write data
    # MATLAB: t = Tiff(out_path, 'r+'); write(t, im); close(t)
    try:
        with rasterio.open(out_path, 'w', **profile) as dst:
            # Ensure the output data type matches the profile
            im_to_write = im.astype(profile['dtype'])
            dst.write(im_to_write, 1)
        print(f"Saved GeoTIFF to {out_path}.")
    except Exception as e:
        warnings.warn(f"Error writing GeoTIFF with rasterio: {e}")

if __name__ == '__main__':
    # Mock usage: requires a 'Cab.tif' file to exist in the exercise/images path
    pass