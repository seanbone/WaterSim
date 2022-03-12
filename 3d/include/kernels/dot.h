#ifndef WATERSIM_DOT_H
#define WATERSIM_DOT_H

#include <immintrin.h>

// returns <x,y> = x.T * y
// Note: expects aligned arrays
/*
double dot(const double *x, const double *y, const unsigned int n) {
	double tmp = 0;
	for(unsigned i = 0 ; i < n; ++i ) {
		tmp += x[i]*y[i];
	}
	return tmp;
}
 */
double dot(const double *x, const double *y, const unsigned int n) {
    double tmp = 0;
    unsigned i = 0;
    // we want 8 FMAs because of Skylake ports lat 4, 2 per cycle
    // so we do 8 doubles per step
    __m256d sol1, sol2, sol3, sol4;
    __m256d vec_x1, vec_x2, vec_y1, vec_y2;
    __m256d vec_x3, vec_x4, vec_y3, vec_y4;
    __m256d vec_x5, vec_x6, vec_y5, vec_y6;
    __m256d vec_x7, vec_x8, vec_y7, vec_y8;
    __m256d vec_por1 = _mm256_setzero_pd();
    __m256d vec_por2 = _mm256_setzero_pd();
    __m256d vec_por3 = _mm256_setzero_pd();
    __m256d vec_por4 = _mm256_setzero_pd();
    __m256d vec_por5 = _mm256_setzero_pd();
    __m256d vec_por6 = _mm256_setzero_pd();
    __m256d vec_por7 = _mm256_setzero_pd();
    __m256d vec_por8 = _mm256_setzero_pd();

    // peel loop for alligned loads for a ! b should be handled too.
    // auto peel = (unsigned long) x & 0x1f;
    // if (peel != 0){
        // peel = ( 32-peel)*d_size_inv;
        // for(; i < peel ;++i){
            // tmp += x[i]*y[i];
        // }
    // }
    //! b should be unalligned loads?
    for(; i < n-32; i +=32 ) {
        vec_x1 = _mm256_load_pd(x+i);
        vec_x2 = _mm256_load_pd(x+i+4);
        vec_x3 = _mm256_load_pd(x+i+8);
        vec_x4 = _mm256_load_pd(x+i+12);
        vec_x5 = _mm256_load_pd(x+i+16);
        vec_x6 = _mm256_load_pd(x+i+20);
        vec_x7 = _mm256_load_pd(x+i+24);
        vec_x8 = _mm256_load_pd(x+i+28);

        vec_y1 = _mm256_load_pd(y+i);
        vec_y2 = _mm256_load_pd(y+i+4);
        vec_y3 = _mm256_load_pd(y+i+8);
        vec_y4 = _mm256_load_pd(y+i+12);
        vec_y5 = _mm256_load_pd(y+i+16);
        vec_y6 = _mm256_load_pd(y+i+20);
        vec_y7 = _mm256_load_pd(y+i+24);
        vec_y8 = _mm256_load_pd(y+i+28);

        // not sure if we get aliasing issues here
        vec_por1 = _mm256_fmadd_pd(vec_x1,vec_y1,vec_por1);
        vec_por2 = _mm256_fmadd_pd(vec_x2,vec_y2,vec_por2);
        vec_por3 = _mm256_fmadd_pd(vec_x3,vec_y3,vec_por3);
        vec_por4 = _mm256_fmadd_pd(vec_x4,vec_y4,vec_por4);
        vec_por5 = _mm256_fmadd_pd(vec_x5,vec_y5,vec_por5);
        vec_por6 = _mm256_fmadd_pd(vec_x6,vec_y6,vec_por6);
        vec_por7 = _mm256_fmadd_pd(vec_x7,vec_y7,vec_por7);
        vec_por8 = _mm256_fmadd_pd(vec_x8,vec_y8,vec_por8);
    }
    // reduce all 32 into 1
    sol1 = _mm256_hadd_pd(vec_por1,vec_por2);
    sol2 = _mm256_hadd_pd(vec_por3,vec_por4);
    sol3 = _mm256_hadd_pd(vec_por5,vec_por6);
    sol4 = _mm256_hadd_pd(vec_por7,vec_por8);

    sol1 = _mm256_hadd_pd(sol1,sol3);
    sol2 = _mm256_hadd_pd(sol2,sol4);
    sol1 = _mm256_hadd_pd(sol1,sol2);
    sol1 = _mm256_hadd_pd(sol1,sol1);
    double sol0[4];
    _mm256_storeu_pd(sol0,sol1);

    tmp += sol0[0] + sol0[2];
    for(; i < n; ++i ) {
        tmp += x[i]*y[i];
    }
    return tmp;
}

#endif
