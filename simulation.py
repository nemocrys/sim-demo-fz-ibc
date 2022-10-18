from objectgmsh import Model, Shape, MeshControlExponential
import gmsh
import os
import numpy as np
import pyelmer.elmer as elmer
from pyelmer.execute import run_elmer_grid, run_elmer_solver
from pyelmer.post import scan_logfile

occ = gmsh.model.occ



####################
# parameters
sim_dir = "./simdata"
current = 100  # A
frequency = 1e6  # Hz
mesh_size_factor = 1  # increase for coarser, decrease for finer mesh
visualize = False  # must be false in docker container

if not os.path.exists(sim_dir):
    os.makedirs(sim_dir)
####################
# geometry modeling
model = Model()

coil_body = occ.add_cylinder(0, 0, 0, 0, 0.008, 0, 0.1)
supply_1 = occ.add_cylinder(0.004 + 0.001, 0.004, 0, 0, 0, 0.2, 0.004)
supply_2 = occ.add_cylinder(-0.004 - 0.001, 0.004, 0, 0, 0, 0.2, 0.004)
coil = occ.fuse([(3, coil_body)], [(3, supply_1), (3, supply_2)])[0][0][1]

slit = occ.add_box(-0.001, 0, 0, 0.002, 1, 1)
hole1 = occ.add_cylinder(0, 0, 0, 0, 0.008, 0, 0.008)
hole2 = occ.add_cone(0, 0.001, 0, 0, 0.007, 0, 0.008, 0.04)
occ.cut([(3, coil)], [(3, slit), (3, hole1), (3, hole2)])

inductor = Shape(model, 3, "inductor", [coil_body])
inductor.mesh_size = 0.0025

air = occ.add_cylinder(0, -0.2, 0, 0, 0.4, 0, 0.2)
air = Shape(model, 3, "air", [air])
air.mesh_size = 0.05

outside = occ.add_cylinder(0, -0.2, 0, 0, 0.4, 0, 1)  # to cut end of power supplies
occ.cut([(3, outside)], air.dimtags, removeTool=False)
occ.cut(inductor.dimtags, [(3, outside)])

occ.cut(air.dimtags, inductor.dimtags, removeTool=False)
air.set_interface(inductor)
model.synchronize()

surf_inductor = Shape(model, 2, "surf_inductor", inductor.get_interface(air))
bnd_air = Shape(model, 2, "bnd_air", [x for x in air.boundaries if x not in air.get_interface(inductor)])
inductor_ends = [x for x in inductor.boundaries if x not in inductor.get_interface(air)]
bnd_supply_1 = Shape(model, 2, "bnd_supply_1", [inductor_ends[0]])
bnd_supply_2 = Shape(model, 2, "bnd_supply_2", [inductor_ends[1]])

model.make_physical()

model.deactivate_characteristic_length()
model.set_const_mesh_sizes()
MeshControlExponential(model, inductor, inductor.mesh_size, exp=1.6, fact=2)

model.generate_mesh(3, optimize="Netgen", size_factor=mesh_size_factor)
if visualize:
    model.show()

model.write_msh(f"{sim_dir}/case.msh")
model.close_gmsh()


####################
# simulation setup
sim = elmer.load_simulation("3D_steady", "./config_elmer.yml")
sim.settings.update({"Angular Frequency": 2*np.pi*frequency})

solver_mgdyn = elmer.load_solver("MGDynamics", sim, "./config_elmer.yml")
solver_mgdyn.data.update({"Angular Frequency": 2*np.pi*frequency})
solver_calc = elmer.load_solver("MGDynamicsCalc", sim, "./config_elmer.yml")
solver_mgdyn.data.update({"Angular Frequency": 2*np.pi*frequency})
solver_output = elmer.load_solver("ResultOutputSolver", sim, "./config_elmer.yml")

eqn_main = elmer.Equation(sim, "eqn_main", [solver_mgdyn, solver_calc], {"name": "eqn_main"})

mat_copper = elmer.load_material("copper-inductor", sim, "./config_elmer.yml")
mat_copper.data.update({"name": "copper-inductor"})
mat_air = elmer.load_material("air", sim, "./config_elmer.yml")
mat_air.data.update({"name": "air"})

inductor = elmer.Body(sim, "inductor", [inductor.ph_id], {"name": "inductor"})
inductor.material = mat_copper
inductor.equation = eqn_main
air = elmer.Body(sim, "air", [air.ph_id], {"name": "air"})
air.material = mat_air
air.equation = eqn_main

bnd_air = elmer.Boundary(sim, "bnd_air", [bnd_air.ph_id], {"name": "bnd_air"})
bnd_air.data.update({
    "AV re {e}": "Real 0.0",
    "AV im {e}": "Real 0.0",
})
bnd_supply_1 = elmer.Boundary(sim, "bnd_supply_1", [bnd_supply_1.ph_id], {"name": "bnd_supply_1"})
bnd_supply_1.data.update({
    "AV re {e}": "Real 0.0",
    "AV im {e}": "Real 0.0",
})
bnd_supply_2 = elmer.Boundary(sim, "bnd_supply_2", [bnd_supply_2.ph_id], {"name": "bnd_supply_2"})
bnd_supply_2.data.update({
    "AV re {e}": "Real 0.0",
    "AV im {e}": "Real 0.0",
    "AV re": "Real 0.0",
    "AV im": "Real 0.0",
})
surf_inductor = elmer.Boundary(sim, "surf_inductor", [surf_inductor.ph_id], {"name": "surf_inductor"})
surf_inductor.data.update({
    "Layer Electric Conductivity": "Real 58.1e+6",
    "Layer Relative Permeability": "Real 1",
})
bnd_current_supply = elmer.Boundary(sim, "bnd_current_supply", data={"name": "bnd_current_supply"})
bnd_current_supply.data.update({
    "Intersection BC(2)": "Integer 1 2",  # TODO parameterize, this is prone for errors!!!
    "Electric Current Density": f"-distribute {current}",
})

####################
# execute simulation
sim.write_sif(sim_dir)
print("Starting ElmerGrid...")
run_elmer_grid(sim_dir, "case.msh")
print("Starting ElmerSolver...")
run_elmer_solver(sim_dir)
err, warn, stats = scan_logfile(sim_dir)
print("Errors:", err)
print("Warnings:", warn)
print("Statistics:", stats)
