#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>

int main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cout << "Usage: ./viewmesh <path/to/meshfile.obj>" << std::endl;
		return 0;
	}

	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	// Load a mesh in OBJ format
	igl::readOBJ(argv[1], V, F);

	// Render the mesh
	igl::opengl::glfw::Viewer viewer;

	viewer.data().set_mesh(V, F);

	viewer.launch();
}
