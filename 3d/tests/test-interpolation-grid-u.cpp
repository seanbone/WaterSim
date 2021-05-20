/*
 * A test to check the grid_interpolate routine of Mac3d.
 */
#include "includes/watersim-test-common.h"
#include "Mac3d.h"

int main() {
	unsigned nx = 5;
	unsigned ny = 5;
	unsigned nz = 5;
	double sx = 10;
	double sy = 10;
	double sz = 10;
	Mac3d mac(nx, ny, nz, sx, sy, sz);
	double* g = mac.pu_;
	double* h = mac.pu_star_;
	double offset_x = -0.5 * mac.cell_sizex_;
	double offset_y = 0;
	double offset_z = 0;

	for (unsigned k = 0; k < nz; k++) {
		for (unsigned j = 0; j < ny; j++) {
			for (unsigned i = 0; i < nx+1; i++) {
				int idx = i + (nx+1)*j + (nx+1)*ny*k;
				mac.pu_     [i + (nx+1) * j + (nx+1) * ny * k] = idx*idx;
				mac.pu_star_[i + (nx+1) * j + (nx+1) * ny * k] = idx*idx;
			}
		}
	}

	// Actual tests

	/******************** LEFT FACE ********************/
	{
		// Linear on top-left edge
		double norm_px = -1;
		double norm_py = 5.5;
		double norm_pz = 1.5;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = 4986; // (54**2 + 84**2)/2

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}

	{
		// Linear on bottom-left edge
		double norm_px = -1;
		double norm_py = -1;
		double norm_pz = 1.5;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = 2250;

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}

	{
		// Linear on front-left edge
		double norm_px = -1;
		double norm_py = 1.5;
		double norm_pz = -1.5;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = 90; // (6**2 + 12**2)/2

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}

	{
		// Linear on back-left edge
		double norm_px = -1;
		double norm_py = 1.5;
		double norm_pz = 5.5;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = 16650; // (126**2 + 132**2)/2

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}

	{
		// Bilinear on left face
		double norm_px = -0.5;
		double norm_py = 0.2;
		double norm_pz = 0.7;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = Mac3d::bilinear_interpolation_normalized(
				g[0], g[30], g[6], g[36], norm_py, norm_pz);

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}




	/******************** RIGHT FACE ********************/
	{
		// Linear on top-right edge
		double norm_px = 5.5;
		double norm_py = 5.5;
		double norm_pz = 3.5;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = 18181; // (119**2 + 149**2)/2

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}

	{
		// Linear on bottom-right edge
		double norm_px = 5.5;
		double norm_py = -0.5;
		double norm_pz = 3.5;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = 12325; // (95**2 + 125**2)/2

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}

	{
		// Linear on front-right edge
		double norm_px = 5.1;
		double norm_py = 1.5;
		double norm_pz = -1;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = 205; // (11**2 + 17**2)/2

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}

	{
		// Linear on back-right edge
		double norm_px = 5.1;
		double norm_py = 1.5;
		double norm_pz = 5.1;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = 17965; // (131**2 + 137**2)/2

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}

	{
		// Bilinear on right face
		double norm_px = 5.5;
		double norm_py = 0.2;
		double norm_pz = 0.7;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = Mac3d::bilinear_interpolation_normalized(
				g[5], g[35], g[11], g[41], norm_py, norm_pz);

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}




	/******************** TOP FACE ********************/

	{
		// Linear on top-front edge
		double norm_px = 1.5;
		double norm_py = 5.1;
		double norm_pz = -1;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = 650.5; // (25**2 + 26**2)/2

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}

	{
		// Linear on top-back edge
		double norm_px = 1.5;
		double norm_py = 5.1;
		double norm_pz = 5.1;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = 21170.5; // (145**2 + 146**2)/2

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}

	{
		// Bilinear on top face
		double norm_px = 1.7;
		double norm_py = 5.2;
		double norm_pz = 1.2;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = Mac3d::bilinear_interpolation_normalized(
				g[55], g[56], g[85], g[86], norm_pz-int(norm_pz), norm_px-int(norm_px));

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}




	/******************** BOTTOM FACE ********************/

	{
		// Linear on front-bottom edge
		double norm_px = 1.5;
		double norm_py = -1;
		double norm_pz = -1;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = 2.5;

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}

	{
		// Linear on back-bottom edge
		double norm_px = 1.5;
		double norm_py = -1;
		double norm_pz = 5.5;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = 14762.5; // (121**2 + 122**2)/2

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}

	{
		// Bilinear on bottom face
		double norm_px = 1.7;
		double norm_py = -0.2;
		double norm_pz = 1.2;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = Mac3d::bilinear_interpolation_normalized(
				g[31], g[32], g[61], g[62], norm_pz-int(norm_pz), norm_px-int(norm_px));

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}




	/******************** FRONT FACE ********************/

	{
		// Bilinear on front face
		double norm_px = 1.5;
		double norm_py = 1.2;
		double norm_pz = -1.7;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = Mac3d::bilinear_interpolation_normalized(
				g[7], g[13], g[8], g[14], norm_px-1, norm_py-1);

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}




	/******************** BACK FACE ********************/

	{
		// Bilinear on back face
		double norm_px = 1.5;
		double norm_py = 1.2;
		double norm_pz = 5.7;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		double expected = Mac3d::bilinear_interpolation_normalized(
				g[127], g[133], g[128], g[134], norm_px-1, norm_py-1);

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}




	/******************** INTERIOR DOMAIN ********************/

	{
		// Inner domain
		double norm_px = 3.76;
		double norm_py = 0.95;
		double norm_pz = 1.25;
		double px = norm_px * mac.cell_sizex_ + offset_x;
		double py = norm_py * mac.cell_sizey_ + offset_y;
		double pz = norm_pz * mac.cell_sizez_ + offset_z;

		unsigned cell_x = norm_px;
		unsigned cell_y = norm_py;
		unsigned cell_z = norm_pz;

		double* p = mac.pu_;

		unsigned i000 = (cell_x    ) + (nx+1) * (cell_y    ) + ny * (nx+1) * (cell_z    );
		unsigned i001 = (cell_x    ) + (nx+1) * (cell_y    ) + ny * (nx+1) * (cell_z + 1);
		unsigned i010 = (cell_x    ) + (nx+1) * (cell_y + 1) + ny * (nx+1) * (cell_z    );
		unsigned i011 = (cell_x    ) + (nx+1) * (cell_y + 1) + ny * (nx+1) * (cell_z + 1);
		unsigned i100 = (cell_x + 1) + (nx+1) * (cell_y    ) + ny * (nx+1) * (cell_z    );
		unsigned i101 = (cell_x + 1) + (nx+1) * (cell_y    ) + ny * (nx+1) * (cell_z + 1);
		unsigned i110 = (cell_x + 1) + (nx+1) * (cell_y + 1) + ny * (nx+1) * (cell_z    );
		unsigned i111 = (cell_x + 1) + (nx+1) * (cell_y + 1) + ny * (nx+1) * (cell_z + 1);

		double expected = Mac3d::trilinear_interpolation_normalized(p[i000], p[i001], p[i010], p[i011],
															        p[i100], p[i101], p[i110], p[i111],
															        norm_px - cell_x, norm_py - cell_y, norm_pz - cell_z);

		double returned = mac.grid_interpolate<Mac3d::GRID_U>(px, py, pz);
		assert(std::abs(expected - returned) < interpolation_tolerance);

		double ret1, ret2;
		std::tie(ret1, ret2) = mac.grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(px, py, pz);
		assert(std::abs(expected - ret1) < interpolation_tolerance);
		assert(std::abs(expected - ret2) < interpolation_tolerance);
	}
}
