#include <igl/readOBJ.h>

#include <igl/opengl/glfw/Viewer.h>



Eigen::MatrixXd V;

Eigen::MatrixXi F;



int main(int argc, char *argv[])

{
      // Load a mesh in OBJ format
      
      igl::readOBJ("../out_meshes/mesh_000300.obj", V, F);
      
      // Plot the mesh
      
      igl::opengl::glfw::Viewer viewer;
      
      viewer.data().set_mesh(V, F);
      
      viewer.launch(); 
}
