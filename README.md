# meteor_pbs18

TODO: merge back into master for submission

Meteorite crashing into the sea - Physically-Based Simulation HS18 project.

Authors: Silvia Nauer, Mikael Stellio, Sean Bone.

# Cloning with submodules

This Git repository has a submodule for libigl. To clone it correctly use either of the following:

    git clone --recurse-submodules git@gitlab.ethz.ch:bones/meteor_pbs18.git
    git clone --recurse-submodules https://gitlab.ethz.ch/bones/meteor_pbs18.git

# Usage

The main simulation is in the folder `3d`. It can be compiled with Make and CMake:

    mkdir -p 3d/build/
    cd 3d/build && cmake ..
    make -j8

This builds two executables, `WaterSim` (the actual simulation) and `ViewMesh` (which can be used to preview a generated mesh).

# 3D Simulation

## Simulation parameters

`WaterSim` has several options. Note that the simulation must be reset before any changes are applied.

 - Display grid: whether to display the MAC grid
 - Export meshes: whether to generate `*.obj` files for each frame (see below)
 - Randomize particles: if checked, the initial particle positions are jittered slightly instead of being on a strict 2x2x2 grid in their starting cell.
 - Max particles display: when running a large simulation, rendering all the particles causes the libigl GUI to lag significantly. Therefore a maximum number of particles actually visualized can be selected. This does not affect the number of particles actually used in simulation.
 - Alpha: the mixing parameter for FLIP and PIC methods. `alpha = 1` means pure PIC, `alpha = 0` means pure FLIP.
 - Timestep: timestep used for simulation, in seconds.
 - The density of the fluid, in kg/m^3
 - Acceleration of gravity in m/s^2
 - Grid resolution XYZ: the number of cells along the corresponding direction (axis)
 - Size XYZ: the size in metres of the simulation environment.

**Caution:** the program expects the grid cells to be cubic in shape, and this assumption is made across the program. So special care sould be taken when setting the simulation size (`sx, sy, sz`) and grid resolution (`nx, ny, nz`) such that `sx/nx = sy/ny = sz/nz`.

## Initial fluid layout

The initial layout of the fluid is determined in `main.cpp`. The function `select_fluid_cells` is used to determine which cells should be initialized as fluid at the beginning of the simulation.
For example:

```
std::vector<bool> select_fluid_cells(size_t nx, size_t ny, size_t nz) {
	std::vector<bool> is_fluid(nx*ny*nz, false);

	for (unsigned k = 0; k < nz; k++) {
		for (unsigned j = 0; j < ny/4; j++) {
			for (unsigned i = 0; i < nx; i++) {
				is_fluid[i + j*nx + nx*ny*k] = true;
			}
		}
    }
    return is_fluid;
}
```

This will fill the lower quarter of the simulation box with fluid.

## Meteor splash

TODO: move splash parameters out of `FLIP.cpp` and write explanation

## Exporting meshes

At the end of each simulation step, if the `Export meshes` option is checked, the level-set function will be evaluated and `igl::copyleft::marching_cubes` is used to generate a mesh for the current frame.
The details of the process can be modified in `3d/include/MeshExporter.h`:

 - `r`: the "radius" of particles, normalized to the cell size (`r = 1` means one cell size)
 - `h`: the radius over which particles are averaged around cell centers, normalized to the cell size (`d = 1` means one cell size)
 - `folder_`: the target folder where meshes will be saved
 - `file_prefix_`: the name to use for files (files are named `<file_prefix_>000000.obj` etc)

**Note on portability**: `MeshExporter::MeshExporter` uses the function [`mkdir`](http://pubs.opengroup.org/onlinepubs/009695399/functions/mkdir.html) to ensure the target directory exists. Should this cause issues on non-POSIX systems, the call can be safely removed from `MeshExporter.cpp` as long as the folder exists.

## Viewing meshes

The exported meshes can of course be imported to other 3D software (Blender, Maya), or it can be previewed with the executable `ViewMesh`. This will open a single file with a simple libigl window. The source for this file is `3d/viewmesh.cpp`.


# 2D simulation

TODO: clean up 2D and write explanation

## Making a video from PNGs
This command will take files numbered `0000.png, 0001.png, 0002.png, ...` in folder `PNG_out` (usually generated in the build directory by WaterSim when "Export PNGs" is enabled) and make an MP4 video of them:
    
    ffmpeg -r 40 -f image2 -s 1280x800 -i PNG_out/%04d.png -vcodec libx264 -crf 20  -pix_fmt yuv420p test.mp4

 - `-r`: frame rate ( = 1/dt, dt is the timestep used in simulation)
 - `-s`: dimensions
 - `-crf`: quality, lower is better
 - `-i`: input images. Change folder name/location if necessary.
