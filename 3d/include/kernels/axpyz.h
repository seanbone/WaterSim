#ifndef WATERSIM_AXPYZ_H
#define WATERSIM_AXPYZ_H

#include <immintrin.h>

//z <- a * x + y
// Note: expects aligned arrays
/*
void axpyz(const unsigned int n, const double a, const double *x, const double *y, double *z) {
    for(unsigned i = 0 ; i < n; ++i ) {
        z[i] = x[i] * a + y[i];
    }
}
 */
void axpyz(const unsigned int n, const double a, const double *x, const double *y, double *z) {
    unsigned i = 0;
    // we want 8 FMAs because of Skylake ports

    __m256d vec_x1, vec_x2;
    __m256d vec_x3, vec_x4;
    __m256d vec_x5, vec_x6;
    __m256d vec_x7, vec_x8;
    __m256d vec_y1, vec_y2;
    __m256d vec_y3, vec_y4;
    __m256d vec_y5, vec_y6;
    __m256d vec_y7, vec_y8;
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
        // peel = ( 32-peel) *d_size_inv;
        // for(; i < peel ;++i){
            // z[i] = x[i] * a + y[i];
        // }
    // }
    for(; i < n-32; i +=32 ) {
        vec_x1 = _mm256_load_pd(x + i);
        vec_x2 = _mm256_load_pd(x + i + 4);
        vec_x3 = _mm256_load_pd(x + i + 8);
        vec_x4 = _mm256_load_pd(x + i + 12);
        vec_x5 = _mm256_load_pd(x + i + 16);
        vec_x6 = _mm256_load_pd(x + i + 20);
        vec_x7 = _mm256_load_pd(x + i + 24);
        vec_x8 = _mm256_load_pd(x + i + 28);
        vec_y1 = _mm256_load_pd(y + i);
        vec_y2 = _mm256_load_pd(y + i + 4);
        vec_y3 = _mm256_load_pd(y + i + 8);
        vec_y4 = _mm256_load_pd(y + i + 12);
        vec_y5 = _mm256_load_pd(y + i + 16);
        vec_y6 = _mm256_load_pd(y + i + 20);
        vec_y7 = _mm256_load_pd(y + i + 24);
        vec_y8 = _mm256_load_pd(y + i + 28);

        // not sure if we get aliasing issues here
        vec_res1 = _mm256_fmadd_pd(vec_a1, vec_x1, vec_y1);
        vec_res2 = _mm256_fmadd_pd(vec_a2, vec_x2, vec_y2);
        vec_res3 = _mm256_fmadd_pd(vec_a3, vec_x3, vec_y3);
        vec_res4 = _mm256_fmadd_pd(vec_a4, vec_x4, vec_y4);
        vec_res5 = _mm256_fmadd_pd(vec_a5, vec_x5, vec_y5);
        vec_res6 = _mm256_fmadd_pd(vec_a6, vec_x6, vec_y6);
        vec_res7 = _mm256_fmadd_pd(vec_a7, vec_x7, vec_y7);
        vec_res8 = _mm256_fmadd_pd(vec_a8, vec_x8, vec_y8);


        _mm256_storeu_pd(z + i, vec_res1);
        _mm256_storeu_pd(z + i + 4, vec_res2);
        _mm256_storeu_pd(z + i + 8, vec_res3);
        _mm256_storeu_pd(z + i + 12, vec_res4);
        _mm256_storeu_pd(z + i + 16, vec_res5);
        _mm256_storeu_pd(z + i + 20, vec_res6);
        _mm256_storeu_pd(z + i + 24, vec_res7);
        _mm256_storeu_pd(z + i + 28, vec_res8);
    }
    for(; i < n; ++i ) {
        z[i] = x[i] * a + y[i];
    }
}


#endif
