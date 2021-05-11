#include "Particles.h"

#include <cstdlib>

Particles::Particles(Particles::particleIdx_t nParticles, const SimConfig &cfg, const Mac3d &macGrid)
		: rcell_size_x_(1.0 / macGrid.get_cell_sizex()), rcell_size_y_(1.0 / macGrid.get_cell_sizey()),
		  rcell_size_z_(1.0 / macGrid.get_cell_sizez()), num_particles_(nParticles) {

	x = (double *) calloc(nParticles, sizeof(double));
	y = (double *) calloc(nParticles, sizeof(double));
	z = (double *) calloc(nParticles, sizeof(double));
	u = (double *) calloc(nParticles, sizeof(double));
	v = (double *) calloc(nParticles, sizeof(double));
	w = (double *) calloc(nParticles, sizeof(double));
}

Particles::~Particles() {
	free(x);
	free(y);
	free(z);
	free(u);
	free(v);
	free(w);
}

void Particles::get_cell_index(particleIdx_t particleIdx,
                                      Mac3d::cellIdx_t &cellIdxX,
                                      Mac3d::cellIdx_t &cellIdxY,
                                      Mac3d::cellIdx_t &cellIdxZ) const {
	cellIdxX = (Mac3d::cellIdx_t) (x[particleIdx] * rcell_size_x_ + 0.5);
	cellIdxY = (Mac3d::cellIdx_t) (y[particleIdx] * rcell_size_y_ + 0.5);
	cellIdxZ = (Mac3d::cellIdx_t) (z[particleIdx] * rcell_size_z_ + 0.5);
}

Particles::particleIdx_t Particles::get_num_particles() const {
	return num_particles_;
}