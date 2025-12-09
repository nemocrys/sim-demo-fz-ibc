# demo-fz-ibc

**Zenodo reference**

Mangetic field and thermal simulation of simple Floating Zone configuration with one-turn flat inductor.

The project is developed and maintained by the [**Model experiments group**](https://www.ikz-berlin.de/en/research/materials-science/section-fundamental-description#c486) at the Leibniz Institute for Crystal Growth (IKZ).

### Referencing

If you use this code in your research, please cite our open-access publication:

> I. Tsiapkinis, A. Wintzer, D. Kaspars, Validation of 3D and 2D thermal and electromagnetic models for high-frequency induction heating in crystal growth processes *Journal of Crystal Growth*,  643 (2024) 127800. [https://doi.org/10.1016/j.jcrysgro.2024.127800](https://doi.org/10.1016/j.jcrysgro.2024.127800).

## Overview

3D numerical simulation for high-frequency induction heating in the FZ process for comparison with magnetic and thermal measurements. The geometry of the experimntal setup is shown below. The current in the experiment and simulation was set to 130A, the frequency is 672kHz.

The simulation uses ElmerFEM with an impedance boundary condition on the conducting surfaces. The heat transfer simulation includes indcution heating, conduction and surface-to-ambient radiation and convection.
Further settings can be found in the configuration file *config_elmer.ylm*.

<p align="center">
  <img src="https://github.com/user-attachments/assets/fe744b13-a79e-4427-828d-f33c5d1d4167" alt="Experimental setup" width="48%">
  <img src="https://github.com/user-attachments/assets/0369c727-7eec-4466-a097-c9381b9e65a3" alt="Setup sketch" width="48%">
</p>

### Configuration


The configuration of the simulation is stored in yml-files. Process parameters (frequency, current) and mesh configuration can be found at the top of *simulation.py*.

The files *config_geo.yml* (geometry configuration), *config_elmer.yml* (simulation setup and material parameters) are usually left unchanged, to make changes to these parameters the section config_update in *config.yml* is used.


### Geometry and simulation setup

The geometry and simulation setup are defined in *simulation.py*.

### Execution

If the required python packages and elmer solver are installed simulations can be executed with the script *simulation.py.

Usage of docker is highly recommended, see next section.

## Computational setup (Docker)

The setup for the simulations is provided in form of a docker image, so just an installation of [Docker](https://docs.docker.com/get-docker/) is required on your system. The image nemocrys/opencgs:v0.3.2 is used (see [opencgs](https://github.com/nemocrys/opencgs) for more information), future versions may require changes to the code.

### Executing the simulations
On Windows, the simulations can be executed with:
```
docker run --rm -v ${PWD}:/home/workdir nemocrys/opencgs:v0.3.2 python3 simulation.py
```

On Linux, the simulations can be executed with:
```
docker run -it --rm -v $PWD:/home/workdir -e LOCAL_UID=$(id -u $USER) -e LOCAL_GID=$(id -g $USER) nemocrys/opencgs:v0.3.2 python3 simulation.py
```

### Running the docker container
On windows the docker container can be started with:
```
docker run --rm -v ${PWD}:/home/workdir nemocrys/opencgs:v0.3.2 bash
```

On Linux, the simulations can be executed with:
```
docker run -it --rm -v $PWD:/home/workdir -e LOCAL_UID=$(id -u $USER) -e LOCAL_GID=$(id -g $USER) nemocrys/opencgs:v0.3.2 bash
```

The simulation can then be run by executing *python3 simulation.py* inside the docker container.

## Results
A selection of simulation results are shown here. For the complete verification and comparison with experiments refer to TODO: Add ref

### Domain and mesh
Sketch of 3D-model simulation domain and (b) detailed view of mesh used in the 3D-model. The mesh consists of approximately 1.5 million elements. The smallest element size is 0.25 mm on the surface of the inductor hole
<p align="center">
  <img src="https://github.com/user-attachments/assets/0a0d8edf-2b0f-4e75-9ef8-47c14306d899" alt="Elmer domain sketch" width="48%">
  <img src="https://github.com/user-attachments/assets/fc4eb37b-5d45-4db0-85de-97a1caaaf12e" width="48%">
</p>

## Results

 Magnetic field lines of 𝐁𝐼𝑚 in the air domain around the inductor, colored with |𝐁| using a logarithmic scale. The surface current distribution |𝐣| on the rod and inductor is also shown. Half of the inductor, cut through the main slit, is shown. Magnetic flux density magnitude distribution

![Magnetic field lines](figures/elmer-results-field-lines.png)

## Power and temperatures
The power is calculated with XX. The total induced power is 8.85 W.
Teh figure below shows the omparison of temperature profiles on the rod surface. The measured profile was extracted from the thermal image. The profile for the 3D-model was taken at 𝛼 = 180◦.
<img width="28%" alt="sim_temperature_htc9" src="https://github.com/user-attachments/assets/2516c6d4-9c62-4659-8660-5abb431869b1" />

## Acknowledgements

[This project](https://nemocrys.github.io/) has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 851768).

<img src="https://github.com/nemocrys/test-cz-induction/blob/main/EU-ERC.png">
