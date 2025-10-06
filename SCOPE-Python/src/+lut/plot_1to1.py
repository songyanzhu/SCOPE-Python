import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import linregress

def plot_1to1(meas, mod, what, out_path, var_name='GPP'):
    """
    // src/+lut/plot_1to1.m
    Generates a 1-to-1 scatter plot with RMSE, Bias, and R^2_adj.
    """
    
    plt.figure()

    # --- 1. Units and Labels ---
    units = r'[W m$^{-2}$]'
    if var_name in ['Actot', 'GPP']:
        units = r'[$\mu$mol CO$_2$ m$^{-2}$ s$^{-1}$]'
        
    xlabel = f"SCOPE {var_name} {units}"
    ylabel = f"{what} {var_name} {units}"

    # --- 2. Scatter Plot ---
    plt.scatter(meas, mod, label='scatter')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # --- 3. 1-to-1 Line (refline(1, 0)) ---
    # Find max/min range for a proper 1-to-1 line
    min_val = np.nanmin([meas, mod])
    max_val = np.nanmax([meas, mod])
    
    plt.plot([min_val, max_val], [min_val, max_val], color='r', linestyle='-', label='1-to-1 line')
    
    # --- 4. Regression Line (refline / fitlm) ---
    # MATLAB: lm = fitlm(meas, mod); 
    # Use scipy for linear regression
    slope, intercept, r_value, p_value, std_err = linregress(meas, mod)
    r_squared = r_value**2
    # R^2_adj is harder to calculate without N and number of predictors (p=1 here)
    n = len(meas)
    p = 1
    r_adj = 1 - (1 - r_squared) * (n - 1) / (n - p - 1)

    # Plot regression line
    plt.plot(meas, intercept + slope * meas, color='g', linestyle='--', label='regression')
    
    # --- 5. Metrics ---
    # MATLAB: rmse = sqrt(mean((meas - mod) .^ 2));
    rmse = np.sqrt(np.mean((meas - mod) ** 2))
    bias = np.mean(mod - meas)
    
    # --- 6. Title and Legend ---
    plt.title(f"{var_name}, RMSE = {rmse:.2f}, bias={bias:.2f}, $R^2_{{adj}}$ = {r_adj:.2f}")
    plt.legend(loc='lower right')
    
    # Ensure axes limits are equal for 1-to-1 comparison
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlim([min_val, max_val])
    plt.ylim([min_val, max_val])
    
    # MATLAB: saveas(gcf, out_path)
    plt.savefig(out_path)
    plt.close()

if __name__ == '__main__':
    # Mock usage example
    meas_mock = np.arange(1, 101) + np.random.randn(100) * 5
    mod_mock = meas_mock * 1.1 + 5 + np.random.randn(100) * 5
    plot_1to1(meas_mock, mod_mock, 'GPR', 'mock_plot_1to1.png', 'Actot')