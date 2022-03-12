#ifndef WATERSIM_AXPYMAXZ_H
#define WATERSIM_AXPYMAXZ_H

#include <immintrin.h>

// z <- a * x + y; returns max |z[i]|
// Note: expects aligned arrays
/*
double axpyzmax(const unsigned int n, const double a, const double *x, const double *y, double *z) {
   double max_abs_val = 0;
   for(unsigned i = 0 ; i < n; ++i ) {
       z[i] = x[i] * a + y[i];
       const double abs_val = std::abs(z[i]);
       if (max_abs_val < abs_val) max_abs_val = abs_val;
   }
   return max_abs_val;
}
*/
double axpyzmax(const unsigned int n, const double a, const double *x, const double *y, double *z) {
    double max_abs_val = 0;
    unsigned i = 0;
    // we want 8 FMAs because of Skylake ports
    // so we do 8 doubles per step
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
    const __m256d vec_a1 = _mm256_set1_pd(a);
    const __m256d vec_a2 = _mm256_set1_pd(a);
    const __m256d vec_a3 = _mm256_set1_pd(a);
    const __m256d vec_a4 = _mm256_set1_pd(a);
    const __m256d vec_a5 = _mm256_set1_pd(a);
    const __m256d vec_a6 = _mm256_set1_pd(a);
    const __m256d vec_a7 = _mm256_set1_pd(a);
    const __m256d vec_a8 = _mm256_set1_pd(a);
    __m256d vec_max1 = _mm256_setzero_pd();
    __m256d vec_max2 = _mm256_setzero_pd();
    __m256d vec_max3 = _mm256_setzero_pd();
    __m256d vec_max4 = _mm256_setzero_pd();
    __m256d vec_max5 = _mm256_setzero_pd();
    __m256d vec_max6 = _mm256_setzero_pd();
    __m256d vec_max7 = _mm256_setzero_pd();
    __m256d vec_max8 = _mm256_setzero_pd();

    // peel loop for alligned loads for a ! b should be handled too.
    // auto peel = (unsigned long) x & 0x1f;
    // if (peel != 0) {
        // peel = (32 - peel) * d_size_inv;
        // for (; i < peel; ++i) {
            // z[i] = y[i] + x[i] * a;
            // const double square_val = z[i] * z[i];
            // if (max_abs_val < square_val) max_abs_val = square_val;
        // }
    // }
    for (; i < n - 32; i += 32) {
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

        vec_res1 = _mm256_mul_pd(vec_res1, vec_res1);
        vec_res2 = _mm256_mul_pd(vec_res2, vec_res2);
        vec_res3 = _mm256_mul_pd(vec_res3, vec_res3);
        vec_res4 = _mm256_mul_pd(vec_res4, vec_res4);
        vec_res5 = _mm256_mul_pd(vec_res5, vec_res5);
        vec_res6 = _mm256_mul_pd(vec_res6, vec_res6);
        vec_res7 = _mm256_mul_pd(vec_res7, vec_res7);
        vec_res8 = _mm256_mul_pd(vec_res8, vec_res8);


        vec_max1 = _mm256_max_pd(vec_res1, vec_max1);
        vec_max2 = _mm256_max_pd(vec_res2, vec_max2);
        vec_max3 = _mm256_max_pd(vec_res3, vec_max3);
        vec_max4 = _mm256_max_pd(vec_res4, vec_max4);
        vec_max5 = _mm256_max_pd(vec_res5, vec_max5);
        vec_max6 = _mm256_max_pd(vec_res6, vec_max6);
        vec_max7 = _mm256_max_pd(vec_res7, vec_max7);
        vec_max8 = _mm256_max_pd(vec_res8, vec_max8);

    }
    for (; i < n; ++i) {
        z[i] = y[i] + x[i] * a;
        const double square_val = z[i] * z[i];
        if (max_abs_val < square_val) max_abs_val = square_val;
    }
    // now find the largest value
    vec_max1 = _mm256_max_pd(vec_max1, vec_max5);
    vec_max2 = _mm256_max_pd(vec_max2, vec_max6);
    vec_max3 = _mm256_max_pd(vec_max3, vec_max7);
    vec_max4 = _mm256_max_pd(vec_max4, vec_max8);

    vec_max1 = _mm256_max_pd(vec_max1, vec_max3);
    vec_max2 = _mm256_max_pd(vec_max2, vec_max4);

    vec_max1 = _mm256_max_pd(vec_max1, vec_max2);

    // reuse max2 to compare same array vals
    // Might be faster to store and compare the doubles directly
    vec_max2 = _mm256_permute_pd(vec_max1, 0xb1);
    vec_max1 = _mm256_max_pd(vec_max2, vec_max1);
    // sanity check: here element 0 and 1 aswell as 2 and 3 should be the same

    // maybe add allignement
    double sol0[4];
    _mm256_storeu_pd(sol0, vec_max1);

    //potential optimisation :
    //const __m256i mask = {1,0,1,0}; // maybe 256
    //_mm256_maskstore_pd(sol0,0xff00ff00,vec_max1);

    return std::sqrt(std::max(std::max(sol0[0], sol0[2]), max_abs_val));
}



#endif
