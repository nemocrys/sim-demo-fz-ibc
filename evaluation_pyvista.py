import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# def cylindrical_to_cartesian(r, phi, z):
#     # COORDINATE TRANSFORMATION SIMULATION -> EXPERIMENT
#     x = r * np.sin(phi)
#     y = z
#     z = r * np.cos(phi)
#     return x, y, z


def readAll_csv_data(filename):
    data = pd.read_csv(filename)
    return data

def virtual_probe(pvd_file, point, probe_bounds = [-5e-3, 5e-3, -5e-3, 5e-3, -5e-3, 5e-3]):
    # box_bounds: size of probe [xmin, xmax, ymin, ymax, zmin, zmax]
    # x, y, z = cylindrical_to_cartesian(r, phi, z)
    box_bounds = [point[0] + probe_bounds[0],
                  point[0] + probe_bounds[1],
                  point[1] + probe_bounds[2],
                  point[1] + probe_bounds[3],
                  point[2] + probe_bounds[4],
                  point[2] + probe_bounds[5]]
    print(box_bounds)
    # Load the PVD file
    
    # pvd_data = pv.read(pvd_file)
    reader = pv.read(pvd_file)
    reader.set_active_time_point(-1)  # last time step (here: steady state iteration)
    mesh = reader.read()
    mesh = mesh.combine()  # combine different blocks (different vtu files) to one mesh
    # # Access the "B" data from the PVD dataset
    # b_data = pvd_data.point_arrays["magnetic flux density re e"]

    # Create a box region of interest (ROI) using the given bounds
    roi = pv.Box(bounds=box_bounds)

    # Extract data within the ROI
    extracted_data = roi.apply(dataset=reader)

    # Calculate the mean value of the extracted "B" data
    mean_value = np.mean(extracted_data.point_arrays["magnetic flux density im e"])

    return mean_value

def evaluate_circle(pvd_file, r, z, resolution=1000, plot=False):
    # create circle
    phi = np.linspace(-np.pi/2, 3/2*np.pi, resolution)
    points = np.array([r*np.sin(phi), z*np.ones(phi.shape), r*np.cos(phi)]).T

    # read pvd
    reader = pv.get_reader(pvd_file)
    reader.set_active_time_point(-1)  # last time step (here: steady state iteration)
    mesh = reader.read()
    mesh = mesh.combine()  # combine different blocks (different vtu files) to one mesh

    # mesh.plot(style="wireframe", scalars="temperature", cmap="turbo")

    # evaluate pvd
    point_set = pv.PointSet(points)
    probe = point_set.sample(point_set)
    values = ((probe.point_data["magnetic flux density re e"]**2 + probe.point_data["magnetic flux density im e"]**2))**0.5
    # values = probe.point_data["magnetic flux density im e"]
    # print(points.shape[0])
    # for i in range(0, points.shape[0]):
    #     virtualProbe = virtual_probe(pvd_file, points[i,:])
    #     print(f'point = {points[i,:]} \t values = {values[i,:]} \t virtualProbe = {virtualprobe.point_data[i,:]}')

    if plot:
        fig, ax = plt.subplots()
        ax.plot(phi, values)
        plt.show()
    
    return phi, values

def evaluate_vertical_line(pvd_file, r, phi, z_start, z_end, resolution=100, plot=False):
    # create line
    z = np.linspace(z_start, z_end, resolution)
    points = np.array([r*np.sin(phi)*np.ones(z.shape), z, r*np.cos(phi)*np.ones(z.shape)]).T

    # read pvd
    reader = pv.get_reader(pvd_file)
    reader.set_active_time_point(-1)  # last time step (here: steady state iteration)
    mesh = reader.read()
    mesh = mesh.combine()  # combine different blocks (different vtu files) to one mesh

    # evaluate pvd
    point_set = pv.PointSet(points)
    probe = point_set.sample(mesh)
    values = probe.point_data["magnetic flux density im e"][:, 1]

    if plot:
        fig, ax = plt.subplots()
        ax.plot(z, values)
        plt.show()
    
    return z, values

def evaluate_line(pvd_file, points, plot=False):
    # read pvd
    reader = pv.get_reader(pvd_file)
    reader.set_active_time_point(-1)  # last time step (here: steady state iteration)
    mesh = reader.read()
    mesh = mesh.combine()  # combine different blocks (different vtu files) to one mesh

    # evaluate pvd
    point_set = pv.PointSet(points)
    probe = point_set.sample(mesh)
    values = ((probe.point_data["magnetic flux density re e"]**2 + probe.point_data["magnetic flux density im e"]**2))**0.5

    if plot:
        fig, ax = plt.subplots()
        ax.plot(points[:,1], values)
        plt.show()
    
    return points, values

def evaluate_crystal_temperature(pvd_file, points, plot=False):
    # read pvd
    reader = pv.get_reader(pvd_file)
    reader.set_active_time_point(-1)  # last time step (here: steady state iteration)
    mesh = reader.read()
    mesh = mesh.combine()  # combine different blocks (different vtu files) to one mesh

    # evaluate pvd
    point_set = pv.PointSet(points)
    probe = point_set.sample(mesh)
    values = probe.point_data['temperature']

    return points, values

def evaluate_circle_temperature(pvd_file, r, z, resolution=1000, plot=False):
    # create circle
    phi = np.linspace(-np.pi/2, 3/2*np.pi, resolution)
    points = np.array([r*np.sin(phi), z*np.ones(phi.shape), r*np.cos(phi)]).T

    # read pvd
    reader = pv.get_reader(pvd_file)
    reader.set_active_time_point(-1)  # last time step (here: steady state iteration)
    mesh = reader.read()
    mesh = mesh.combine()  # combine different blocks (different vtu files) to one mesh

    # mesh.plot(style="wireframe", scalars="temperature", cmap="turbo")

    # evaluate pvd
    point_set = pv.PointSet(points)
    probe = point_set.sample(mesh)
    values = probe.point_data["temperature"]
    # values = ((probe.point_data["magnetic flux density re e"]**2 + probe.point_data["magnetic flux density im e"]**2))**0.5

    if plot:
        fig, ax = plt.subplots()
        ax.plot(phi, values)
        plt.show()
    
    return phi, values

def return_sum_nodal_heat(pvd_file, combine=False):
    # read pvd
    reader = pv.get_reader(pvd_file)
    mesh = reader.read()
    if combine:
        mesh = mesh.combine()  # combine different blocks (different vtu files) to one mesh
    # print(mesh.array_names)
    field = mesh["nodal joule heating"]      # replace with your field name
    sum = field.sum()

    cell_mesh = mesh.point_data_to_cell_data()
    field_cell = cell_mesh["nodal joule heating"]
    volumes = cell_mesh.compute_cell_sizes()["Volume"]
    integral = np.sum(field_cell * volumes)
    # print("Integral:", integral)
    # print("Raw sum:", sum)

    return sum

def compute_heat_flux_from_vtu(vtu_file, k):
    """
    Load a VTU file, convert point->cell data,
    compute temperature gradient and boundary heat flux.

    Parameters
    ----------
    vtu_file : str
        Path to .vtu file.
    k : float
        Thermal conductivity [W/(m·K)].

    Returns
    -------
    total_flux : float
        Total outward heat flux [W].
    q_surf : np.ndarray
        Heat flux vector at boundary cells (shape: n_boundary_cells x 3).
    """

    # --- Load file ---
    mesh = pv.read(vtu_file)

    # --- convert point data to cell data ---
    cell_mesh = mesh.point_data_to_cell_data()

    # make sure the temperature field exists
    if "temperature" not in cell_mesh.cell_data:
        raise KeyError(
            f"Temperature field 'temperature' not found in VTU. "
            f"Available cell arrays: {list(cell_mesh.cell_data.keys())}"
        )

    # --- compute gradient of temperature (cell data) ---
    # use compute_derivative (vtkGradientFilter wrapper)
    grad_mesh = cell_mesh.compute_derivative(
        scalars="temperature",
        gradient="gradT",    # name of resulting gradient array
        preference="cell",   # we want the cell-data field
    )

    # gradient is now a cell array (n_cells x 3)
    gradT = grad_mesh.cell_data["gradT"]

    # --- compute heat flux vector: q = -k * gradT ---
    q = -k * gradT
    grad_mesh.cell_data["q"] = q

    # --- extract external surface ---
    surf = grad_mesh.extract_surface()

    # compute cell normals on the surface
    surf = surf.compute_normals(cell_normals=True, point_normals=False)

    # get normals and flux as *cell data* on the surface
    normals = surf.cell_data["Normals"]   # shape: (n_surf_cells, 3)
    q_surf = surf.cell_data["q"]          # shape: (n_surf_cells, 3)

    # --- normal heat flux qn = q · n ---
    qn = np.einsum("ij,ij->i", q_surf, normals)

    # --- get each face area ---
    areas = surf.compute_cell_sizes()["Area"]   # length n_surf_cells

    # --- integrate heat flux (outward) ---
    total_flux = np.sum(qn * areas)

    return total_flux, q_surf, qn

if __name__ == "__main__":
    pvd_file = "simdata_size-factor=2_discontinuous-bodies/results/case.pvd"

    y = 0.04  # height  (relative to vessel bottom)
    r = 0.055  # radius
    phi = 0
    resolution = 100
    # evlauate_circle(pvd_file, r, y, resolution, True)

    evaluate_vertical_line(pvd_file, 0, 0, 0.05, 0.15, 100, True)
