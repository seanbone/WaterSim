#ifndef WATERSIM_XMAX_H
#define WATERSIM_XMAX_H

#include <immintrin.h>

// returns max |x[i]| for i in 0:n-1
// Note: expects aligned arrays
/*
double xmax(const int n, const double* x) {
   double max_abs_val = 0;
   for (int i = 0; i < n; i++) {
       const double abs_val = std::abs(x[i]);
       if (max_abs_val < abs_val) max_abs_val = abs_val;
   }
   return max_abs_val;
}
 */
double xmax(const int n, const double* x) {
    double max_abs_val = 0;
    unsigned i = 0;
    // we want 8 Mults because of Skylake ports
    __m256d vec_x1, vec_x2;
    __m256d vec_x3, vec_x4;
    __m256d vec_x5, vec_x6;
    __m256d vec_x7, vec_x8;
    __m256d vec_max1 = _mm256_setzero_pd();
    __m256d vec_max2 = _mm256_setzero_pd();
    __m256d vec_max3 = _mm256_setzero_pd();
    __m256d vec_max4 = _mm256_setzero_pd();
    __m256d vec_max5 = _mm256_setzero_pd();
    __m256d vec_max6 = _mm256_setzero_pd();
    __m256d vec_max7 = _mm256_setzero_pd();
    __m256d vec_max8 = _mm256_setzero_pd();

    // peel loop for alligned loads for x
    // unsigned peel = (unsigned long) x & 0x1f;
    // if (peel != 0) {
        // peel = (32 - peel) * d_size_inv;
        // for (; i < peel; ++i) {
            // const double square_val = x[i]*x[i];
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

        vec_x1 = _mm256_mul_pd(vec_x1, vec_x1);
        vec_x2 = _mm256_mul_pd(vec_x2, vec_x2);
        vec_x3 = _mm256_mul_pd(vec_x3, vec_x3);
        vec_x4 = _mm256_mul_pd(vec_x4, vec_x4);
        vec_x5 = _mm256_mul_pd(vec_x5, vec_x5);
        vec_x6 = _mm256_mul_pd(vec_x6, vec_x6);
        vec_x7 = _mm256_mul_pd(vec_x7, vec_x7);
        vec_x8 = _mm256_mul_pd(vec_x8, vec_x8);

        vec_max1 = _mm256_max_pd(vec_x1,vec_max1);
        vec_max2 = _mm256_max_pd(vec_x2,vec_max2);
        vec_max3 = _mm256_max_pd(vec_x3,vec_max3);
        vec_max4 = _mm256_max_pd(vec_x4,vec_max4);
        vec_max5 = _mm256_max_pd(vec_x5,vec_max5);
        vec_max6 = _mm256_max_pd(vec_x6,vec_max6);
        vec_max7 = _mm256_max_pd(vec_x7,vec_max7);
        vec_max8 = _mm256_max_pd(vec_x8,vec_max8);
    }
    for (; i < n; ++i) {
        const double square_val = x[i]*x[i];
        if (max_abs_val < square_val) max_abs_val = square_val;
    }
    // now find the largest value
    vec_max1 = _mm256_max_pd(vec_max1,vec_max5);
    vec_max2 = _mm256_max_pd(vec_max2,vec_max6);
    vec_max3 = _mm256_max_pd(vec_max3,vec_max7);
    vec_max4 = _mm256_max_pd(vec_max4,vec_max8);

    vec_max1 = _mm256_max_pd(vec_max1,vec_max3);
    vec_max2 = _mm256_max_pd(vec_max2,vec_max4);

    vec_max1 = _mm256_max_pd(vec_max1,vec_max2);

    // reuse max2 to compare same array vals
    // Might be faster to store and compare the doubles directly
    vec_max2 = _mm256_permute_pd(vec_max1, 0xb1);
    vec_max1 = _mm256_max_pd(vec_max2, vec_max1);

    double sol0 [4];
    _mm256_storeu_pd(sol0,vec_max1);
    return std::sqrt(std::max(std::max(sol0[0], sol0[2]), max_abs_val));
}


#endif
