import numpy as np
import os
from PIL import Image

# Mock external functions from the same package
# In a real setup, these would be imported from other .py files in src/lut
from .plot_image import plot_image 

# Assumed external function (if needed for utility)
def imread(path):
    """Mock for imread - reads an image into a numpy array"""
    if path.endswith('.tif'):
        try:
            import rasterio
            with rasterio.open(path) as src:
                return src.read(1)
        except ImportError:
            # Fallback for mock environment
            return np.ones((100, 100)) 
    else:
        return np.ones((100, 100))

def change_detection(rmse_im_path=None, gpr_im_path=None):
    """
    // src/+lut/change_detection.m
    """
    if rmse_im_path is None:
        in_dir = os.path.join('..', 'exercise')
        rmse_im_path = os.path.join(in_dir, 'results_rmse.tif')
        gpr_im_path = os.path.join(in_dir, 'results_gpr.tif')

    out_path = os.path.join(os.path.dirname(rmse_im_path), 'results_rmse_gpr.png')
    
    assert os.path.exists(rmse_im_path), \
        f"Did not find `{rmse_im_path}` image.\nHave you retrieved on full_set.csv with lut.use_rmse()?"
    assert os.path.exists(gpr_im_path), \
        f"Did not find `{gpr_im_path}` image.\nHave you retrieved on full_set.csv with lut.use_gpr()?"
    
    # imread function needs proper raster/gdal handling for GeoTIFFs
    lut_im = imread(rmse_im_path)
    gpr_im = imread(gpr_im_path)
    
    im_dif = lut_im - gpr_im
    
    # lut.plot_image(im_dif, '(RMSE - GPR)', out_path)
    plot_image(im_dif, '(RMSE - GPR)', out_path)

if __name__ == '__main__':
    change_detection()