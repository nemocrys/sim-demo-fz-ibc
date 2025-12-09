import numpy as np
import pandas as pd

def read_exp_data(filename, target_r, target_z, tolerance=1e-3):
    data = pd.read_csv(filename, delim_whitespace=True)

    data = data.applymap(lambda x: np.nan if x > 1e37 else x)

    filtered_data = data[(abs(data["r[mm]"] - target_r) <= tolerance) & (abs(data["z[mm]"] - target_z) <= tolerance)]
    return filtered_data

def read_exp_data2(filename, target_alpha, target_z, tolerance=1e-3):
    data = pd.read_csv(filename, delim_whitespace=True)

    data = data.applymap(lambda x: np.nan if x > 1e37 else x)

    filtered_data = data[(abs(data["a[deg]"] - target_alpha) <= tolerance) & (abs(data["z[mm]"] - target_z) <= tolerance)]
    return filtered_data

def readAll_exp_data(filename):
    data = pd.read_csv(filename, delim_whitespace=True)

    data = data.applymap(lambda x: np.nan if x > 1e37 else x)

    # filtered_data = data[abs(data["r[mm]"] - target_r) <= tolerance]
    return data

def find_unique_r_values(filename):
    data = pd.read_csv(filename, delim_whitespace=True)
    unique_r_values = data["r[mm]"].unique()
    return unique_r_values.tolist()

def find_unique_z_values(filename):
    data = pd.read_csv(filename, delim_whitespace=True)
    unique_z_values = data["z[mm]"].unique()
    return unique_z_values.tolist()

if __name__ == "__main__":
    exp_file = "expdata/results_z12.txt"
    rvals = find_unique_r_values(exp_file)
    data = read_exp_data(exp_file, rvals[4])
    print(data)
