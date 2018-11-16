# PBS project: meteor

Simulation parameters: FPS, # frames, grid size xyz, grid resolution xyz, solid configuration, initial fluid configuration

Program pipeline:
1. Construct MAC grid
    - size & resolution
    - 3 arrays for velocities $u$, $v$, $w$
    - one array for pressures
    - one bool array to store solid cells
    - precompute $A$ pressure matrix diagonal
2. Spawn initial FLIP particles
	- 2x2x2 subgrid of each cell
3. Perform FLIP for duration of 1 frame
	- Particle-based advection for intermediate velocities: $\frac{\partial \vec u}{\partial t} +\vec u \cdot \nabla\vec u = \vec f$
		- Forward Euler to update particle positions
		- Particle to grid tran
	- Solve for pressures: $\nabla\cdot\vec u = 0$ and $\nabla^2p = \nabla\cdot\vec u$
	- Compute final velocities $\vec u^{n+1}$: $\frac{\partial \vec u}{\partial t} = -\nabla p$
4. Compute level set function
5. Run marching cubes to generate surface mesh
6. Repeat from 3. until all frames have been generated
7. Import mesh sequence (frames) into Maya & render


# FLIP algorithm
In each iteration:
1. Forward Euler (/RK2) to update particle positions
2. Transfer particle velocities to grid $\rightarrow u^*$ & save copy
3. Apply external forces (gravity etc) with Forward Euler: $u^* += \vec g\Delta t$
4. Construct pressure matrix $A$
5. Compute right-hand side (divergence) $d$
6. Solve $Ap = d$ with MICCG(0)
7. Update grid velocities with pressure gradients $\rightarrow u^{n+1}$
8. Update particle velocities by mixing FLIP and PIC
9. Rinse&repeat
