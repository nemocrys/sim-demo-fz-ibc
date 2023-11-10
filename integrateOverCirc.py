import csv
import math
import matplotlib.pyplot as plt
import numpy as np
def calculate_distance(x, y):
    return math.sqrt(x**2 + y**2)

datafile = "simdata_130A_665kHz_slice.csv"  # Replace with the actual filename
points = []
with open(datafile, "r") as file:
    reader = csv.DictReader(file)
    max_distance = 0
    max_point = None

    for row in reader:
        point_id = int(row["Point ID"])
        x = float(row["Points_0"])
        y = float(row["Points_1"])
        distance = calculate_distance(x+0.005, y-0.004)
        # print(distance)
        if distance>0.0038:
            points.append(np.array([x+0.005,y-0.004,float(row["surface current im_2"]),float(row["surface current re_2"])]))
        if distance > max_distance:
            max_distance = distance
            max_point = (point_id, x, y)

    print("Outer-most Point:")
    print("Point ID:", max_point[0])
    print("Coordinates:", max_point[1], max_point[2])

points = np.array(points)
print(points)


x = points[:, 0]
y = points[:, 1]
value1 = points[:, 2]
value2 = points[:, 3]

# Sort the data points based on the x-coordinate
distances = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
sort_indices = np.argsort(x)
x = x[sort_indices]
y = y[sort_indices]
print(np.arctan(x/y))
value1 = value1[sort_indices]
value2 = value2[sort_indices]
fig,ax =plt.subplots(1,3)
ax[0].scatter(x,y,c=value1)
ax[0].axis('equal')

# Compute the length of the line
line_length = np.sqrt(np.diff(x)**2 + np.diff(y)**2).sum()

# Integrate value1 and value2 using the trapezoidal rule
integral_value1 = np.trapz(value1, x) / line_length
integral_value2 = np.trapz(value2, x) / line_length

print("Integral of value1:", integral_value1)
print("Integral of value2:", integral_value2)

distances = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
print(np.sqrt(x**2 + y**2))
print([x,y])
cumulative_distances = np.cumsum(distances)
cumulative_distances = np.insert(cumulative_distances, 0, 0)  # Add the starting point

# Normalize the cumulative distances to the range [0, 2π]
normalized_distances = cumulative_distances * (2*np.pi / cumulative_distances[-1])

# Calculate the variable s
s = np.arctan2(x,y)* 0.004
sort_indices = np.argsort(s)
s = s[sort_indices]
x = x[sort_indices]
y = y[sort_indices]
print(np.arctan(x/y))
value1 = value1[sort_indices]
value2 = value2[sort_indices]
integral_value1 = np.trapz(value1, s) 
integral_value2 = np.trapz(value2, s) 
integral_value3 = np.trapz(np.sqrt(value2**2+value1**2), s) 
# print("Variable s:", s)
print("Integral of value1:", integral_value1)
print("Integral of value2:", integral_value2)
print("Integral of abs:", integral_value3)
# ax[1].plot(s, value1)
ax[1].plot(s/0.004, value1)
ax[1].plot(s/0.004, value2)
ax[1].set_xlabel("phi")
ax[1].set_ylabel("surface current density")
ax[1].legend(["surface current density Im","surface current density Re"])
plt.tight_layout()

fig,ax =plt.subplots()
ax.plot(s/0.004+np.pi, value1)
ax.plot(s/0.004+np.pi, value2)
ax.set_xlabel("phi")
ax.set_ylabel("surface current density")
ax.legend(["surface current density Im","surface current density Re"])
plt.tight_layout()
print(92*np.sqrt(2))
plt.show()