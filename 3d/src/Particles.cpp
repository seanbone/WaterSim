#include "Particles.h"

#include <cstdlib>

Particles::Particles(Particles::particleIdx_t nParticles, const Mac3d &macGrid)
		: rcell_size_x_(1.0 / macGrid.get_cell_sizex()), rcell_size_y_(1.0 / macGrid.get_cell_sizey()),
		  rcell_size_z_(1.0 / macGrid.get_cell_sizez()), num_particles_(nParticles) {

	x =  (double *) calloc(nParticles, sizeof(double));
	y =  (double *) calloc(nParticles, sizeof(double));
	z =  (double *) calloc(nParticles, sizeof(double));
	u =  (double *) calloc(nParticles, sizeof(double));
	v =  (double *) calloc(nParticles, sizeof(double));
	w =  (double *) calloc(nParticles, sizeof(double));

    /*
    x =  new (std::align_val_t(32)) double [nParticles];
    y =  new (std::align_val_t(32)) double [nParticles];
    z =  new (std::align_val_t(32)) double [nParticles];
    u =  new (std::align_val_t(32)) double [nParticles];
    v =  new (std::align_val_t(32)) double [nParticles];
    w =  new (std::align_val_t(32)) double [nParticles];
    */
}

Particles::~Particles() {
	free(x);
	free(y);
	free(z);
	free(u);
	free(v);
	free(w);
}