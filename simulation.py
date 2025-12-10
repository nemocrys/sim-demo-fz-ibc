from objectgmsh import Model, Shape, MeshControlExponential
import gmsh
import os
import numpy as np
import pyelmer.elmer as elmer
import yaml
from pyelmer.execute import run_elmer_grid, run_elmer_solver
from pyelmer.post import scan_logfile

occ = gmsh.model.occ


####################
# parameters

with open("config_geo.yml") as f:
        config_geo = yaml.safe_load(f)

with open("config_elmer.yml") as f:
        config_elmer = yaml.safe_load(f)

sim_dir = "./simdata"

## sim inputs
current = 1.838/20*1000*np.sqrt(2)  #
frequency = 672000  # Hz
heat_transfer_coefficient = 9.58 # W/m^2/K
mesh_size_factor = 2  # global increase for coarser, decrease for finer mesh

## ibc conductivities
layer_cond_copper = 58.1e6
layer_cond_tin = 4.38e6

visualize = False  # must be false in docker container

feed_shift = 1.e-3 # moves feed to the top, if positive

if not os.path.exists(sim_dir):
    os.makedirs(sim_dir)

####################
# mesh sizes
inductor_mesh_size = 0.001
mesh_inductor_refinement = 4
feed_mesh_size = 0.001
air_mesh_size = 0.01

inductor_mesh_exp = 1.6
inductor_mesh_fact = 2

feed_mesh_exp = 1.6
feed_mesh_fact = 2

####################
# Inductor geometry, single turn flat inductor with power supplies
r_inner = config_geo["fz_inductor"]["r_inner"]
r_outer = config_geo["fz_inductor"]["r_outer"]
r_outerCoil = config_geo["fz_inductor"]["r_outerCoil"]
r_bevel = config_geo["fz_inductor"]["r_bevel"]
h_inner = config_geo["fz_inductor"]["h_inner"]
h_outer = config_geo["fz_inductor"]["h_outer"]
l_supply = config_geo["fz_inductor"]["l_supply"]
t_slit = config_geo["fz_inductor"]["t_slit"]
r_coil = h_outer/2
r_conus = np.polyfit(np.array([h_inner, h_outer]),np.array([r_inner, r_bevel]),1)[1]

# dummy rod geometry
X0_feed = config_geo["fz_crystal"]["X0"]
r_feed = config_geo["fz_crystal"]["r"]
l_feed = config_geo["fz_crystal"]["l"]

X0_air = config_geo["surrounding"]["X0"]
r_air = config_geo["surrounding"]["r"]
h_air = config_geo["surrounding"]["h"]

####################
# geometry modeling
model = Model()

coil_body = occ.add_cylinder(0, 0, 0, 0, h_outer, 0, r_outer)
supply_1 = occ.add_cylinder(r_coil + t_slit/2, r_coil, 0, 0, 0, l_supply, r_coil)
supply_2 = occ.add_cylinder(-r_coil - t_slit/2, r_coil, 0, 0, 0, l_supply, r_coil)
coil_torus = occ.add_torus(0, r_coil, 0, r_outer, r_coil, zAxis = [0, 1, 0])

# model.synchronize()
# model.show()

coil = occ.fuse([(3, coil_body)], [(3, supply_1), (3, supply_2), (3, coil_torus)])[0][0][1]

slit = occ.add_box(-t_slit/2, 0, 0, t_slit, 1, 1)
hole1 = occ.add_cylinder(0, 0, 0, 0, h_outer, 0, r_inner)
hole2 = occ.add_cone(0, 0.00, 0, 0, h_outer, 0, r_conus, r_bevel)
occ.cut([(3, coil)], [(3, slit), (3, hole1), (3, hole2)])

inductor = Shape(model, 3, "inductor", [coil_body])
inductor.mesh_size = inductor_mesh_size

feed = occ.add_cylinder(X0_feed[0], X0_feed[1], 0, 0, l_feed, 0, r_feed)
feed = Shape(model, 3, "feed", [feed])
feed.mesh_size = feed_mesh_size

air = occ.add_cylinder(X0_air[0], X0_air[1], 0, 0, h_air, 0, r_air)
air = Shape(model, 3, "air", [air])
air.mesh_size = air_mesh_size

outside = occ.add_cylinder(X0_air[0], X0_air[1], 0, 0, h_air, 0, r_air*5)  # to cut end of power supplies
occ.cut([(3, outside)], air.dimtags, removeTool=False)
occ.cut(inductor.dimtags, [(3, outside)])

occ.cut(air.dimtags, inductor.dimtags + feed.dimtags, removeTool=False)
air.set_interface(inductor)
air.set_interface(feed)

model.synchronize()
# model.show() 

surf_inductor = Shape(model, 2, "surf_inductor", inductor.get_interface(air))
surf_feed = Shape(model, 2, "surf_feed", feed.get_interface(air))
bnd_air = Shape(model, 2, "bnd_air", [x for x in air.boundaries if x not in air.get_interface(inductor) + air.get_interface(feed)])
inductor_ends = [x for x in inductor.boundaries if x not in inductor.get_interface(air)]
bnd_supply_1 = Shape(model, 2, "bnd_supply_1", [inductor_ends[0]])
bnd_supply_2 = Shape(model, 2, "bnd_supply_2", [inductor_ends[1]])

inductor_inner_ring = Shape(model, 2, "inductor_inner_ring", [29])  # helper boundary for mesh refinement, extracted manually
# feed_bot = Shape(model, 2, "feed_bot", [21])  # helper boundary for mesh refinement, extracted manually

model.make_physical()


model.deactivate_characteristic_length()
model.set_const_mesh_sizes()
MeshControlExponential(model, inductor, inductor.mesh_size, exp=inductor_mesh_exp, fact=inductor_mesh_fact)
MeshControlExponential(model, feed, feed.mesh_size, exp=feed_mesh_exp, fact=feed_mesh_fact)
MeshControlExponential(model, inductor_inner_ring, inductor.mesh_size/mesh_inductor_refinement, exp=1.6, fact=2)

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
solver_heat = elmer.load_solver("HeatSolver", sim, "./config_elmer.yml")
solver_output = elmer.load_solver("ResultOutputSolver", sim, "./config_elmer.yml")
solver_post = elmer.load_solver("calculateIntegralValues", sim, "./config_elmer.yml")

eqn_main = elmer.Equation(sim, "eqn_main", [solver_mgdyn, solver_calc, solver_heat], {"name": "eqn_main"})

mat_copper = elmer.load_material("copper-ibc", sim, "./config_elmer.yml")
mat_copper.data.update({"name": "copper-ibc"})
mat_air = elmer.load_material("air", sim, "./config_elmer.yml")
mat_air.data.update({"name": "air"})
mat_tin = elmer.load_material("tin-ibc", sim, "./config_elmer.yml")
mat_tin.data.update({"name": "tin-ibc"})

inductor = elmer.Body(sim, "inductor", [inductor.ph_id], {"name": "inductor", "inductor body": "Logical True"})
inductor.material = mat_copper
inductor.equation = eqn_main
air = elmer.Body(sim, "air", [air.ph_id], {"name": "air", "air body": "Logical True"})
air.material = mat_air
air.equation = eqn_main

feed = elmer.Body(sim, "feed", [feed.ph_id], {"name": "feed", "feed body": "Logical True"})
feed.material = mat_tin
feed.equation = eqn_main

joule_heat = elmer.BodyForce(sim, "joule_heat")
joule_heat.data = {
    "Temperature Load": 'Equals "Nodal Joule Heating"'
}
feed.body_force = joule_heat
inductor.body_force = joule_heat

bnd_air = elmer.Boundary(sim, "bnd_air", [bnd_air.ph_id], {"name": "bnd_air"})
bnd_air.data.update({
    "AV re {e}": "Real 0.0",
    "AV im {e}": "Real 0.0",
    "Temperature": 293.15,
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
    "Layer Electric Conductivity": f"Real {layer_cond_copper}",
    "Layer Relative Permeability": "Real 1",
    "Temperature": 293.15,
    "Save Scalars": "Logical True",
})
surf_feed = elmer.Boundary(sim, "surf_feed", [surf_feed.ph_id], {"name": "surf_feed"})
surf_feed.data.update({
    "Layer Electric Conductivity": f"Real {layer_cond_tin}",
    "Layer Relative Permeability": "Real 1",
    "Heat transfer coefficient": heat_transfer_coefficient,
    "External Temperature": 293.15,
    "Radiation": "Idealized",
    "Save Scalars": "Logical True",
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
