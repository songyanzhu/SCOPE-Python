import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.linear_model import LinearRegression

def plot_output_one2one():
    """
    Reads flux data from two specified SCOPE output files and generates
    a 1:1 plot comparison with fit metrics (RMSE, Bias, R^2 adjusted).
    
    NOTE: The MATLAB code has multiple commented-out path_y assignments. 
    The last one is used.
    """
    
    # Define paths
    # These paths are hardcoded in the MATLAB example and assume a specific directory structure.
    path_x = os.path.join('output', 'SCOPE_sparse_2020-06-07-1346', 'fluxes.csv')
    # path_y is set to the last non-commented value in the MATLAB snippet
    path_y = os.path.join('output', 'bigleaf_my_2020-06-08-0030', 'fluxes.csv')
    
    try:
        # Read dataframes
        df_x = pd.read_csv(path_x)
        df_y = pd.read_csv(path_y)
    except FileNotFoundError as e:
        print(f"Error: One or both files not found. Ensure paths are correct and files exist. {e}")
        return

    # List of flux variables to plot
    flu_names = ['Rnctot', 'lEctot', 'Hctot', 'Actot', 'Tcave', 
                 'Rnstot', 'lEstot', 'Hstot', 'Gtot', 'Tsave',
                 'Rntot', 'lEtot', 'Htot', 'rss']

    # Setup the figure
    fig, axes = plt.subplots(3, 5, figsize=(18, 10))
    axes = axes.flatten() # Flatten the 2D array of axes for easy indexing

    for i, flux in enumerate(flu_names):
        if i >= len(axes):
            break # Safety break if flu_names exceeds subplot count

        ax = axes[i]
        
        # Extract data
        if flux not in df_x.columns or flux not in df_y.columns:
            print(f"Warning: Flux '{flux}' not found in one or both dataframes. Skipping.")
            continue
            
        x = df_x[flux].values
        y = df_y[flux].values
        
        # Find non-NaN indices
        not_nan_mask = ~np.isnan(x) & ~np.isnan(y)
        x_nn = x[not_nan_mask]
        y_nn = y[not_nan_mask]

        # Scatter plot
        ax.plot(x, y, 'o', markerfacecolor='r', markeredgecolor='r', alpha=0.6)
        
        # --- Linear Model and Fit ---
        if len(x_nn) > 1:
            # Linear regression: y = m*x + c
            # NumPy's polyfit is equivalent to the MATLAB `polyfit` for degree 1
            lm = np.polyfit(x_nn, y_nn, 1)
            
            # Predict line span
            predict_x = np.array([np.nanmin(x), np.nanmax(x)])
            fit = np.polyval(lm, predict_x)
            
            # Plot the fitted line
            ax.plot(predict_x, fit, 'r:', linewidth=2.5)

            # --- Metrics ---
            # RMSE
            rmse = np.sqrt(np.nanmean((x - y) ** 2))
            
            # Bias (Mean difference)
            bias = np.nanmean(x - y)
            
            # R^2 Adjusted (Requires explicit linear model fitting, no simple polyfit metric)
            # Using scikit-learn for R-squared adjusted calculation
            X_model = x_nn.reshape(-1, 1)
            model = LinearRegression().fit(X_model, y_nn)
            y_pred = model.predict(X_model)
            
            # Calculate R^2
            ss_res = np.sum((y_nn - y_pred) ** 2)
            ss_tot = np.sum((y_nn - np.mean(y_nn)) ** 2)
            r_squared = 1 - (ss_res / ss_tot)
            
            # R^2 Adjusted (r2_adj = 1 - (1-R^2) * (N-1) / (N-P-1))
            N = len(x_nn) # Number of observations
            P = 1 # Number of predictors (x)
            if N > P + 1:
                r_squared_adj = 1 - (1 - r_squared) * (N - 1) / (N - P - 1)
            else:
                r_squared_adj = np.nan
                
            # Title with metrics (using rich text for red color via LaTeX-like syntax for reference)
            title_text = f"{flux}\n{rmse:.2g} (rmse), {bias:.2g} (bias), R^2adj: {r_squared_adj:.1g}"
            ax.set_title(title_text)
            
        else:
            # Not enough data for fit
            ax.set_title(f"{flux}\n(Insufficient data)")
            
        ax.set_xlabel('SCOPE')
        ax.set_ylabel('bigleaf my')
        
        # 1:1 line (refline(1, 0) in MATLAB)
        min_val = np.nanmin([x, y])
        max_val = np.nanmax([x, y])
        line_range = np.array([min_val, max_val])
        ax.plot(line_range, line_range, 'k--', linewidth=1.0) # Black dashed line

    # Remove unused subplots
    for j in range(len(flu_names), len(axes)):
        fig.delaxes(axes[j])

    # Super title
    fig.suptitle('Bigleaf ebal no Fc vs original SCOPE', fontsize=16)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust layout to make room for suptitle
    plt.show()

# To run this function, you would need the specified files in the correct directories
# plot_output_one2one()