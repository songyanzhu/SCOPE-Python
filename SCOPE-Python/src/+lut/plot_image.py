import numpy as np
import matplotlib.pyplot as plt
import os

def plot_image(im, what, out_path):
    """
    // src/+lut/plot_image.m
    Generates an image plot with a colorbar and title.
    """
    
    plt.figure()
    
    # MATLAB: imagesc(im) is equivalent to imshow with a color map
    # MATLAB: caxis([quantile(im(:), 0.01), quantile(im(:), 0.99)])
    
    # 1. Determine color limits using 1st and 99th percentiles
    im_flat = im.flatten()
    im_flat = im_flat[~np.isnan(im_flat)]
    
    if im_flat.size > 0:
        v_min, v_max = np.nanpercentile(im_flat, [1, 99])
    else:
        v_min, v_max = 0, 1
        
    # 2. Plot the image
    img = plt.imshow(im, cmap='viridis', vmin=v_min, vmax=v_max)
    
    # 3. Colorbar and Titles
    cb = plt.colorbar(img)
    cb.set_label(r'[$\mu$mol CO$_2$ m$^{-2}$ s$^{-1}$]')
    plt.title(f"GPP with {what}")

    # 4. Save
    # MATLAB: saveas(gcf, out_path)
    plt.savefig(out_path)
    plt.close()

if __name__ == '__main__':
    # Mock usage example
    mock_im = np.random.rand(100, 100) * 50
    mock_im[5:10, 5:10] = np.nan
    plot_image(mock_im, 'RMSE', 'mock_plot_image.png')