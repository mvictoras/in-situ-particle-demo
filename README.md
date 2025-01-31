# MPI Ascent Particle Simulation

This repository contains a minimal example of a **particle simulation** using **MPI**, **Conduit**, and **Ascent** for in-situ visualization on HPC systems. Particles are initialized randomly inside a sphere of radius \(R\), move with a fixed velocity magnitude \(u\), and “bounce” off the sphere’s boundary. Each timestep, simulation data is published to Ascent for visualization.

## Features

- **MPI Parallel**: Particles are split across MPI ranks.  
- **Conduit-Based**: Particle data is stored in Conduit [Blueprint](https://llnl-conduit.readthedocs.io/en/latest/blueprint.html) format.  
- **Ascent In-Situ Visualization**: Pseudocolor plots of velocity magnitude are generated and saved as images each timestep.  
- **Sphere Boundary Reflection**: Particles reflect off a spherical boundary to remain within the domain.

## Build Instructions

### Configure and Build
```bash
mkdir build
cd build
cmake -DAscent_DIR=<path to ascent> ..
make
```

### Run
```bash
mpirun -n 2 ./particle_sim <X> <R> <T> <u>
```

Where:
- **X**: Total number of particles (must be divisible by the number of ranks).
- **R**: Sphere radius.
- **T**: Number of timesteps.
- **u**: Particle velocity magnitude.

For example:
```bash
mpirun -n 2 ./particle_sim 10000 10.0 50 0.1
```
This runs on 2 MPI ranks with 10,000 total particles (5,000 per rank), a sphere radius of 10.0, 50 timesteps, and velocity magnitude of 0.1.

## Example Visualization Output
![Example Visualization Output](example.png "Example Visualization Output")


## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
