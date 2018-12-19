# meteor_pbs18

Meteorite crashing into the sea - Physically-Based Simulation HS18 project.

Authors: Silvia Nauer, Mikael Stellio, Sean Bone.

This is a 3D FLIP solver implemented as part of the course Physically-Based Simulation for Computer Graphics (ETHZ autumn semester 2018), the goal of which was to create a video of a meteorite crashing into the sea.

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

To produce the effect of a meteorite hitting the water, we apply radial forces within a fixed distance of a moving point, effectively simulating an object in motion. This is implemented in `FLIP::explode`, which is called from `FLIP::step_FLIP`:

    explode(dt, step, 15, 0, 15, 2, 800);

The first two parameters are the timestep and current step number. The following three arguments are the coordinated (in cells) of the point the "meteorite" is directed at. The last two parameters are the radius (in cells) of the meteorite and the force it applies.

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

The folder `2d` contains the 2D version of our simulation. This was the blueprint for the final 3D version. This means the code is not quite as "clean" and there might be some unfixed bugs still lurking, but it also has a few features that weren't carried over to the 3D version:

 - Show pressure field: checking this will color the background to reflect the pressures calculated by the simulation, dark red being the highest pressure in this time step, and dark blue the lowest. We opted for a per-cell coloring as opposed to a smooth one, as this is closer to what the program was actually calculating and proved to be more useful for debugging.
 - Display velocity arrows: this is another feature we used for debugging. If checked, lines are drawn on the grid representing the U and V velocity fields calculated on the MAC grid. This proved invaluable to debugging particle-to-grid and grid-to-particle transfers, as well as velocity extrapolation.
- Export PNGs: if checked, this will generate a PNG of the particle layout at each step, and save it in `PNG_out`. These can then be composed into a video (see below). The image resolution can be set with "PNG size xy" and a maximum number of images to be generated can be imposed with "Max PNGs". A word of caution: this functionality is implemented in the visualization loop (not the simulation thread), so it may cause the GUI to lag when enabled.

**Note on portability**: `WaterSim::exportPNG` uses the function [`mkdir`](http://pubs.opengroup.org/onlinepubs/009695399/functions/mkdir.html) to ensure the target directory exists. Should this cause issues on non-POSIX systems, the call can be safely removed from `WaterSim.cpp` as long as the folder exists.

**Caution:** similarly to the 3D version, this program assumes the mesh cells to be square, so care should be taken when setting XY resolution and size.

## Making a video from PNGs
This command uses the program `ffmpeg` to take files numbered `0000.png, 0001.png, 0002.png, ...` in folder `PNG_out` (for instance generated by WaterSim when "Export PNGs" is enabled) and make an MP4 video of them:
    
    ffmpeg -r 40 -f image2 -s 1280x800 -i PNG_out/%04d.png -vcodec libx264 -crf 20  -pix_fmt yuv420p myvid.mp4

 - `-r`: frame rate ( = 1/dt, dt is the timestep used in simulation)
 - `-s`: dimensions
 - `-crf`: quality, lower is better
 - `-i`: input images. Change folder name/location if necessary.


# Known issues & bugs

 - We originally intended to implement interaction with static grid-aligned solids. Due to time constraints, we didn't finish implementing them, since they weren't required for our video. However, there is still some code referencing solids: this has no effect on the simulation, and calls to `Mac3d::is_solid` will evaluate to false. Since we may yet come back to this project (it's fun!), we decided to leave in what we had done.
 - Our implementation of the level-set function is still somewhat rudimentary, and as a result the reconstructed meshes are imperfect. This produces some slight instabilities in the surface. This is especially visible in small splashes of just a few particles, where the droplets can be seen "flickering" and changing size between frames.
 - The original video we presented also has some issues. For aethetics, we introduced some small-scale displacements to the water surface in Blender. However, as the meteor first enters the water, this causes the water's surface to wriggle unnaturally. This is probably due to the deformation of the mesh changing the procedural displacements. Unfortunately, we realised this too late to re-render the entire video for the presentation.

# Future features

If/when we get back to this project, there's some things we'd like to implement next.

 - Optimization: even though the compiler helps, there's still plenty of room for improvement. The solver could also be parallelized to and extent to better use computational resources.
 - Smoke (for the meteor)
 - Foam (it looks cool)
 - Some kind of open-boundary conditions, allowing waves and splashes to disappear out of the simulation instead of being forced to remain in the box
 - Particle sources/sinks
 - Grid-aligned stationary solids
