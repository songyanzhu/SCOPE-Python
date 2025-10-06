import pandas as pd
import numpy as np
import os

def pick100(full_set_path=None, n_pixels=100):
    """
    // src/+lut/pick100.m
    Randomly selects a subset of pixels from the full dataset for validation.
    """
    if full_set_path is None:
        full_set_path = os.path.join('..', 'exercise', 'full_set.csv')
        
    out_file = os.path.join(os.path.dirname(full_set_path), 'validation_in.csv')
    
    assert not os.path.exists(out_file), f"`{out_file}` file already exists, delete it first"
    
    # MATLAB: df = readtable(full_set_path);
    df = pd.read_csv(full_set_path)
    
    # MATLAB: ind = randi(size(df, 1), n_pixels, 1);
    n_rows = df.shape[0]
    # np.random.randint(low, high, size): low is inclusive, high is exclusive. 
    # Use np.random.choice for sampling without replacement (better for non-randomness in randi)
    ind_0based = np.random.choice(n_rows, size=n_pixels, replace=False)
    
    # MATLAB: df_out = df(ind, :);
    df_out = df.iloc[ind_0based]
    
    # MATLAB: writetable(df_out, out_file)
    df_out.to_csv(out_file, index=False)
    print(f"saved {n_pixels} pixels to `{out_file}`")

if __name__ == '__main__':
    pick100()