import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import yaml
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib.ticker import FuncFormatter


from evaluation_pyvista import evaluate_crystal_temperature, evaluate_circle_temperature, return_sum_nodal_heat

colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]
font = {'size'   : 12}
matplotlib.rc('font', **font)

#### input
simulation_dir = "."
prefix = "res"

#### calibration factors
per = [2,0,1]
####

### sim data
with open(f"{simulation_dir}/config_geo.yml") as f:
        config_geo = yaml.safe_load(f)
X0_feed = config_geo["fz_crystal"]["X0"]
r_feed = config_geo["fz_crystal"]["r"]
l_feed = config_geo["fz_crystal"]["l"]

pvd_file = f"{simulation_dir}/simdata/results/case.pvd"

z = np.linspace(X0_feed[1]*1, (X0_feed[1]+l_feed)*1, 1000)
points = np.array([0*z, z, 0*z+r_feed*0.99]).T # at supply
points2 = np.array([0*z, z, 0*z-r_feed*0.99]).T # at 180°
points, T = evaluate_crystal_temperature(pvd_file, points)
points2, T2 = evaluate_crystal_temperature(pvd_file, points2)
T = T-273.15
T2 = T2-273.15

fig,ax = plt.subplots(figsize=(3.8,5))
ax.plot(T2, z*1000,color='tab:blue', label="3D-model")

# ax.set_title("T on rod surface")
ax.set_xlabel(f"$T$ [°C]")   
ax.set_ylabel(f"$z$ [mm]")   
ax.set_ylim([-100,100])
ax.grid()
legend = ax.legend(fontsize=12, loc="lower left")
fig.tight_layout()
fig.savefig(f"temperature.png", dpi=600)

fig,ax = plt.subplots(figsize=(3.8,5))

z = np.linspace(X0_feed[1]*1, (X0_feed[1]+l_feed)*1, 1000)
points = np.array([0*z, z, 0*z+r_feed*0.99]).T # at supply
points2 = np.array([0*z, z, 0*z-r_feed*0.99]).T # at 180°
points, T = evaluate_crystal_temperature(pvd_file, points)
points2, T2 = evaluate_crystal_temperature(pvd_file, points2)
T = T-273.15
T2 = T2-273.15

print(T.max())
print(T2.max())
print(T.max()-T2.max())
# window_size = 15
# smoothed_signal = np.convolve(max_values_df.to_numpy().flatten(), np.ones(window_size)/window_size, mode='valid')
# ax.plot(smoothed_signal, z_exp2[7:-7]*M-23,color='k', label="measured")
ax.plot(T, z*1000,'--',color='tab:blue', label="3D-model a=0")
ax.plot(T2, z*1000,color='tab:blue', label="3D-model a=180")

# ax.set_title("T on rod surface")
ax.set_xlabel(f"$T$ [°C]")   
ax.set_ylabel(f"$z$ [mm]")   
ax.set_ylim([-100,100])
ax.grid()
legend = ax.legend(fontsize=12, loc="lower left")
fig.tight_layout()
fig.savefig(f"temperature-2.png", dpi=600)


phi, T_over_alpha = evaluate_circle_temperature(pvd_file, r_feed*0.99, 0, resolution=1000, plot=False)
fig,ax = plt.subplots()
ax.plot(phi*180/np.pi, T_over_alpha-273.15)

ax.set_xlabel("$\\alpha$ [°]")
ax.set_ylabel(f"$T$ [°C]")
ax.grid()
plt.tight_layout()
fig.savefig("temperature-azimuthal.png")

# Elmer:
x_values = T2
y_values = z * 1000
idx45 = np.abs(y_values - 4.5).argmin()
idx72 = np.abs(y_values - 72).argmin()
idx0 = np.abs(y_values - 0).argmin()
idx18 = np.abs(y_values - (-18)).argmin()
idx30 = np.abs(y_values - (-30)).argmin()
print(f"T(z={y_values[idx45]} = {x_values[idx45]}")
print(f"T(z={y_values[idx72]} = {x_values[idx72]}")
print(f"T(z={y_values[idx0]} = {x_values[idx0]}")
print(f"T(z={y_values[idx18]} = {x_values[idx18]}")
print(f"T(z={y_values[idx30]}) = {x_values[idx30]}")
# Find the index of the maximum x value
max_x_index = np.argmax(x_values)
# Retrieve the corresponding x and y values
max_x_value = x_values[max_x_index]
max_y_value = y_values[max_x_index]
# print(f"Maximum x value: {max_x_value}")
# print(f"Corresponding y value: {max_y_value}")
print(f"Tmax(z={max_y_value}) = {max_x_value}")

## radiation-to-ambient:
Tamb = 293.15
epsilon = 0.064
htc = 9.58
kappa = 60
sigmaB = 5.67074419e-8

elmer_T = T2 + 273.15

elmer_heatflux = (sigmaB * epsilon * (elmer_T**4 - Tamb**4) + htc*(elmer_T-Tamb))
elmer_radiative = (sigmaB * epsilon * (elmer_T**4 - Tamb**4))

# print(np.trapz(elmer_heatflux,z))

print(f'Qcond = {np.trapz(elmer_heatflux,z)*2*np.pi*r_feed}')

print(f'Qrad = {np.trapz(elmer_radiative,z)*2*np.pi*r_feed}')

print(f'Sum nodal joule heat feed  = {return_sum_nodal_heat("simdata/results/case_feed_t0002.vtu")} W')
print(f'Sum nodal joule heat inductor  = {return_sum_nodal_heat("simdata/results/case_inductor_t0002.vtu")} W')
print(f'Sum nodal joule heat total  = {return_sum_nodal_heat(pvd_file, True)} W')
