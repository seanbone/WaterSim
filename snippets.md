# Plotting points
See method `renderRenderGeometry` in `demos/1_cannonball/CannonBallSim.h`.
Use an `igl::opengl::glfw::Viewer::data()` to add point coordinates and set point size/colour:

    viewer.data().add_points(
        m_trajectories[trajectory][point].transpose(),
        m_trajectoryColors[trajectory]);
    
    viewer.data().point_size = 5;

