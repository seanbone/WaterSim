# meteor_pbs18

Meteorite crashing into the sea - Physically-Based Simulation HS18 project.

Authors: Silvia Nauer, Mikael Stellio, Sean Bone.

# Cloning with submodules
This Git repository has a submodule for libigl. To clone it correctly use either of the following:

    git clone --recurse-submodules git@gitlab.ethz.ch:bones/meteor_pbs18.git
    git clone --recurse-submodules https://gitlab.ethz.ch/bones/meteor_pbs18.git

# TODO
 - [ ] Set up 2D visualization for FLIP particles & cells
 - [ ] 2D FLIP working correctly
   - [ ] MAC2D data structure
   - [ ] Updates to velocities & pressures
 - [ ] Extend to 3D FLIP && 3D viz
   - [ ] MAC3D data structure
   - [ ] Updates to velocities & pressures
 - [ ] Marching cubes & export mesh at each frame
 - [ ] Import meshes into Maya/Blender for rendering

# Program pipeline
Simulation parameters: FPS, # frames, grid size xyz, grid resolution xyz, solid configuration, initial fluid configuration

Program pipeline:
1. Construct MAC grid
    - size & resolution
    - 3 arrays for velocities $u$, $v$, $w$
    - one array for pressures
    - one bool array to store solid cells
    - precompute $A$ pressure matrix diagonal
3. Spawn initial FLIP particles (?)
4. Perform FLIP for duration of 1 frame
	- Particle-based advection for intermediate velocities: $\frac{\partial \vec u}{\partial t} +\vec u \cdot \nabla\vec u = \vec f$
	- Solve for pressures: $\nabla\cdot\vec u = 0$ and $\nabla^2p = \nabla\cdot\vec u$
	- Compute final velocities $\vec u^{n+1}$: $\frac{\partial \vec u}{\partial t} = -\nabla p$
5. Compute level set function
6. Run marching cubes to generate surface mesh
7. Repeat from 3. until all frames have been generated
8. Import mesh sequence (frames) into Maya & render