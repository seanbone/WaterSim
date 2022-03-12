#ifndef WATERSIM_AXY_H
#define WATERSIM_AXY_H

#include <immintrin.h>

// y <- a * x
// Note: expects aligned arrays
/*
void axy(const unsigned int n, const double a, const double *x, double *y) {
	for(unsigned i = 0 ; i < n; ++i ) {
		y [i] = x[i] * a;
	}
}
 */
void axy(const unsigned int n, const double a, const double *x, double *y) {
    unsigned i = 0;
    // we want Mul has lat 4 and 2 ports on skylake
    // https://www.agner.org/optimize/instruction_tables.pdf p 275.
    // so we do 8 doubles per step
    __m256d vec_x1, vec_x2;
    __m256d vec_x3, vec_x4;
    __m256d vec_x5, vec_x6;
    __m256d vec_x7, vec_x8;
    __m256d vec_res1, vec_res2;
    __m256d vec_res3, vec_res4;
    __m256d vec_res5, vec_res6;
    __m256d vec_res7, vec_res8;
    const __m256d  vec_a1 = _mm256_set1_pd(a);
    const __m256d  vec_a2 = _mm256_set1_pd(a);
    const __m256d  vec_a3 = _mm256_set1_pd(a);
    const __m256d  vec_a4 = _mm256_set1_pd(a);
    const __m256d  vec_a5 = _mm256_set1_pd(a);
    const __m256d  vec_a6 = _mm256_set1_pd(a);
    const __m256d  vec_a7 = _mm256_set1_pd(a);
    const __m256d  vec_a8 = _mm256_set1_pd(a);


    // peel loop for alligned loads for a ! b should be handled too.
    // auto peel = (unsigned long) x & 0x1f;
    // if (peel != 0){
        // peel = ( 32-peel) * d_size_inv;
        // for(; i < peel ;++i){
            // y [i] = x[i] * a;
        // }
    // }
    for(; i < n-32; i +=32 ) {
        vec_x1      = _mm256_load_pd(x+i);
        vec_x2      = _mm256_load_pd(x+i+4);
        vec_x3      = _mm256_load_pd(x+i+8);
        vec_x4      = _mm256_load_pd(x+i+12);
        vec_x5      = _mm256_load_pd(x+i+16);
        vec_x6      = _mm256_load_pd(x+i+20);
        vec_x7      = _mm256_load_pd(x+i+24);
        vec_x8      = _mm256_load_pd(x+i+28);

        vec_res1 =  _mm256_mul_pd(vec_x1,vec_a1);
        vec_res2 =  _mm256_mul_pd(vec_x2,vec_a2);
        vec_res3 =  _mm256_mul_pd(vec_x3,vec_a3);
        vec_res4 =  _mm256_mul_pd(vec_x4,vec_a4);
        vec_res5 =  _mm256_mul_pd(vec_x5,vec_a5);
        vec_res6 =  _mm256_mul_pd(vec_x6,vec_a6);
        vec_res7 =  _mm256_mul_pd(vec_x7,vec_a7);
        vec_res8 =  _mm256_mul_pd(vec_x8,vec_a8);
        _mm256_storeu_pd(y+i,vec_res1);
        _mm256_storeu_pd(y+i+4,vec_res2);
        _mm256_storeu_pd(y+i+8,vec_res3);
        _mm256_storeu_pd(y+i+12,vec_res4);
        _mm256_storeu_pd(y+i+16,vec_res5);
        _mm256_storeu_pd(y+i+20,vec_res6);
        _mm256_storeu_pd(y+i+24,vec_res7);
        _mm256_storeu_pd(y+i+28,vec_res8);
    }
    for(; i < n; ++i ) {
        y [i] = x[i] * a;
    }
}


#endif
