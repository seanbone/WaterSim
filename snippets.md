# Plotting points
See method `renderRenderGeometry` in `demos/1_cannonball/CannonBallSim.h`.
Use an `igl::opengl::glfw::Viewer::data()` to add point coordinates and set point size/colour:

    viewer.data().add_points(
        m_trajectories[trajectory][point].transpose(),
        m_trajectoryColors[trajectory]);
    
    viewer.data().point_size = 5;

# Plotting segments
See `Gui::drawArrow`:

    m_viewer.data_list[0].add_edges(arrow.start, arrow.end, arrow.color);

Use `m_viewer.append_mesh()` to add a new mesh.

# Marching cubes
`libigl` has a marching cubes implementation, see `lib/libigl/include/igl/copyleft/marching_cubes.h`.

# References
- FLIP explanation: [http://www.danenglesson.com/](http://www.danenglesson.com/)
- [Libigl tutorial](https://libigl.github.io/tutorial/)
- 
