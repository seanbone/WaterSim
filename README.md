# flip_asl21: Optimization of a FLIP algorithm

<img height="250" src="https://gitlab.inf.ethz.ch/COURSE-ASL/asl21/team32/-/raw/ad459c11b081bf09a87802cb96c4803e54b58cf1/docs/dam-break-low-res.png"></img>
<img height="250" src="https://gitlab.inf.ethz.ch/COURSE-ASL/asl21/team32/-/raw/ad459c11b081bf09a87802cb96c4803e54b58cf1/docs/dam-break-high-res.png"></img>

Fluid simulation via FLIP (Fluid Implicit Particle) Method.

Physically-Based Simulation HS18 project. Optimized for better performance as part of the Advanced Systems Lab FS21 project.

Original authors: Silvia Nauer, Mikael Stellio, Sean Bone.

Optimizations by: Christoph Amevor, Felix Illes, Mikael Stellio, Sean Bone.

This is a 3D FLIP solver originally implemented as part of the course Physically-Based Simulation for Computer Graphics (ETHZ autumn semester 2018), the goal of which was to create a video of a meteorite crashing into the sea.
The project was revisited with a new team as part of the course Advanced Systems Lab (ETHZ spring semester 2021) with the objective of improving single-core numerical performance.

Images: Visualization of speedup of the optimized version (bottom) over the original implementation (top) by comparison of simulations of different sizes requiring an equal amount of run time.

# Cloning with submodules

This Git repository has a submodule for libigl. To clone it correctly use either of the following:

    git clone --recurse-submodules git@gitlab.inf.ethz.ch:COURSE-ASL/asl21/team32.git
    git clone --recurse-submodules https://gitlab.inf.ethz.ch/COURSE-ASL/asl21/team32.git

For more details on the dependencies for `libigl`, check out the [`libigl` documentation](https://libigl.github.io/tutorial/).
### Note for linux users

Many linux distributions do not include `gcc` and the basic development tools in their default installation. On Ubuntu, you need to install the following packages:

```
sudo apt-get install cmake make build-essential libx11-dev mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev libxrandr-dev libxi-dev libxmu-dev libblas-dev libxinerama-dev libxcursor-dev
```

If you are using linux with a virtual machine on Windows, it is *recommended* to use **Visual Studio** instead.

### Note for Windows users

`libigl` supports the **Microsoft Visual Studio 2015** compiler and later, in *64bit* mode. You can download *Visual Studio 2019 Community* for free from [here](https://visualstudio.microsoft.com/vs/).

# Compiling and running

The main simulation is in the folder `3d`. It can be compiled with Make and CMake:

    mkdir -p 3d/build/
    cd 3d/build && cmake ..
    make -j8

This builds three executables:
* `watersim-gui` starts a GUI with the simulation.
* `watersim-cli` runs without a GUI, and therefore with fewer dependencies. It requires a configuration file as input. Run `./watersim-cli -h` for help.
* `viewmesh`, which can be used to preview an OBJ mesh.

## Reference data

The reference data used in the unit tests for validation of all the sub-steps of 
the FLIP::step_FLIP() function is stored in the `3d/validation_data` and 
comprises of the initial configuration file `validation-config.json` and the 
netCDF file `ref.nc`, where the model state at a specific timestep is stored.

To generate a new `ref.nc` file in the build directory simply enable the CMake 
option `WRITE_REFERENCE` as follows:

    mkdir -p 3d/build/
    cd 3d/build && cmake -DWRITE_REFERENCE=ON ..
    make -j8

Then run either `watersim-gui` or `watersim-cli` at least until the timestep specified in `3d/CMakeLists.txt` is reached.

Note that to avoid writing reference files in further runs the `WRITE_REFERENCE` 
option needs to be explicitly turned off by re-running 
`cmake -DWRITE_REFERENCE=OFF ..`.

## Unit tests

The `3d/tests` directory defines unit tests for the 3D case.
* To run the tests manually, execute `make watersim-tests` in the build directory to build and run tests.
* If you're using an IDE to run the tests, you can use the `watersim-tests-build` target to build all tests without running them, and then use CTest to run them.
* To add a test, simply add a `*.cpp` file to the `3d/tests` folder. CMake will set up a test for each file it finds in that folder (non-recursively).

# 3D Simulation

## Simulation parameters

`watersim` has several options. Note that in GUI mode, the simulation must be reset before any changes to settings are applied (except the first two, which only affect the viewport and not the simulation itself).

 - Show grid: whether to display the MAC grid
 - Max particles display: when running a large simulation, rendering all the particles causes the libigl GUI to lag significantly. Therefore a maximum number of particles actually visualized can be selected. This does not affect the number of particles actually used in simulation.
 - Export meshes: whether to generate `*.obj` files for each frame (see below)
 - Randomize particles: if checked, the initial particle positions are jittered slightly instead of being on a strict 2x2x2 grid in their starting cell.
 - Meteor force: if checked, the `FLIP::explode` method will be called at each step, generating the "meteor splash" effect.
 - Alpha: the mixing parameter for FLIP and PIC methods. `alpha = 1` means pure PIC, `alpha = 0` means pure FLIP.
 - Timestep: timestep used for simulation, in seconds.
 - Max steps: how many steps to run the simulation for. If negative, no limit is set.
 - The density of the fluid, in kg/m^3
 - Acceleration of gravity in m/s^2
 - Grid resolution XYZ: the number of cells along the corresponding direction (axis)
 - Size XYZ: the size in metres of the simulation environment.
 - Fluid region: used to select a region to be filled with fluid at the start of the simulation. It is specified by two points (coordinates in meters) which define an axis-aligned bounding box. All cells whose center lies in this region are flagged as fluid.

**Caution:** the program expects the grid cells to be cubic in shape, and this assumption is made across the program. So special care sould be taken when setting the simulation size (`sx, sy, sz`) and grid resolution (`nx, ny, nz`) such that `sx/nx = sy/ny = sz/nz`.

### Configuration files
Simulation settings can be exported and loaded to/from JSON files.
* In CLI mode, a configuration file is required. By default, it will attempt to read `3d/config.json`, but a different value can be passed as command-line argument. Since `3d/config.json` is ignored by git, you will need to provide this file manually when first cloning the project (or generate it with the GUI).
* In GUI mode, the `3d/config.json` is read. If it is not found, a set of defaults (specified in `SimConfig::setDefaults`) is loaded. Using the GUI buttons, you can save the current settings to `3d/config.json`, reload the settings from `3d/config.json`, or reset to defaults.

Here is an example `config.json` file:
```
{
    "alpha": 0.01,
    "applyMeteorForce": false,
    "density": 1000.0,
    "displayGrid": false,
    "exportMeshes": false,
    "fluidRegion": [
        [
            22.0,
            22.0,
            22.0
        ],
        [
            90.0,
            110.0,
            90.0
        ]
    ],
    "gravity": 9.81,
    "gridResolution": [
        15,
        15,
        15
    ],
    "jitterParticles": true,
    "maxParticlesDisplay": 42420,
    "maxSteps": -1,
    "systemSize": [
        120.0,
        120.0,
        120.0
    ],
    "timeStep": 0.025
}
```

## Exporting meshes

At the end of each simulation step, if the `Export meshes` option is checked, the level-set function will be evaluated and `igl::copyleft::marching_cubes` is used to generate a mesh for the current frame.
The details of the process can be modified in `3d/include/MeshExporter.h`:

 - `r`: the "radius" of particles, normalized to the cell size (`r = 1` means one cell size)
 - `h`: the radius over which particles are averaged around cell centers, normalized to the cell size (`d = 1` means one cell size)
 - `folder_`: the target folder where meshes will be saved
 - `file_prefix_`: the name to use for files (files are named `<file_prefix_>000000.obj` etc)

**Note on portability**: `MeshExporter::MeshExporter` uses the function [`mkdir`](http://pubs.opengroup.org/onlinepubs/009695399/functions/mkdir.html) to ensure the target directory exists. Should this cause issues on non-POSIX systems, the call can be safely removed from `MeshExporter.cpp` as long as the folder exists.

## Viewing meshes

The exported meshes can of course be imported to other 3D software (Blender, Maya), or it can be previewed with the executable `viewmesh`. This will open a single file with a simple libigl window. The source for this file is `3d/viewmesh.cpp`.

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

# Future features

If/when we get back to this project, there's some things we'd like to implement next.
 - Parallelization: Apart from the single-core optimizations, several parts (e.g. the pressure solver) could also be parallelized to better use computational resources.
 - Smoke (for the meteor)
 - Foam (it looks cool)
 - Some kind of open-boundary conditions, allowing waves and splashes to disappear out of the simulation instead of being forced to remain in the box
 - Particle sources/sinks
 - Grid-aligned stationary solids
