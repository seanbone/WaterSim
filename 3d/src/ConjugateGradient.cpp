#include "ConjugateGradient.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

#define square(X)	((X) * (X))

// ********* Kernels **********

// variables for optimisation
const unsigned d_size_inv = 1/sizeof(double);

// returns x.T * y
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
    auto peel = (unsigned long) x & 0x1f;
    if (peel != 0){
        peel = ( 32-peel)*d_size_inv;
        for(; i < peel ;++i){
            tmp += x[i]*y[i];
        }
    }
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

// y <- a * x
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
    auto peel = (unsigned long) x & 0x1f;
    if (peel != 0){
        peel = ( 32-peel) * d_size_inv;
        for(; i < peel ;++i){
            y [i] = x[i] * a;
        }
    }
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

// y <- a * x + y
/*
void axpy(const unsigned int n, const double a, const double *x, double *y) {
	for(unsigned i = 0 ; i < n; ++i ) {
		y [i] += x[i] * a;
	}
}
*/

void axpy(const unsigned int n, const double a, const double *x, double *y) {
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
    const __m256d  vec_a1 = _mm256_set1_pd(a);
    const __m256d  vec_a2 = _mm256_set1_pd(a);
    const __m256d  vec_a3 = _mm256_set1_pd(a);
    const __m256d  vec_a4 = _mm256_set1_pd(a);
    const __m256d  vec_a5 = _mm256_set1_pd(a);
    const __m256d  vec_a6 = _mm256_set1_pd(a);
    const __m256d  vec_a7 = _mm256_set1_pd(a);
    const __m256d  vec_a8 = _mm256_set1_pd(a);

    // peel loop for alligned loads for a ! b should be handled too.
    auto peel = (unsigned long) x & 0x1f;
    if (peel != 0){
        peel = ( 32-peel) *d_size_inv;
        for(; i < peel ;++i){
            y [i] += x[i] * a;
        }
    }
    for(; i < n-32; i +=32 ) {
        vec_x1      = _mm256_load_pd(x+i);
        vec_x2      = _mm256_load_pd(x+i+4);
        vec_x3      = _mm256_load_pd(x+i+8);
        vec_x4      = _mm256_load_pd(x+i+12);
        vec_x5      = _mm256_load_pd(x+i+16);
        vec_x6      = _mm256_load_pd(x+i+20);
        vec_x7      = _mm256_load_pd(x+i+24);
        vec_x8      = _mm256_load_pd(x+i+28);
        vec_y1      = _mm256_load_pd(y+i);
        vec_y2      = _mm256_load_pd(y+i+4);
        vec_y3      = _mm256_load_pd(y+i+8);
        vec_y4      = _mm256_load_pd(y+i+12);
        vec_y5      = _mm256_load_pd(y+i+16);
        vec_y6      = _mm256_load_pd(y+i+20);
        vec_y7      = _mm256_load_pd(y+i+24);
        vec_y8      = _mm256_load_pd(y+i+28);

        // not sure if we get aliasing issues here
        vec_res1 = _mm256_fmadd_pd(vec_a1,vec_x1,vec_y1);
        vec_res2 = _mm256_fmadd_pd(vec_a2,vec_x2,vec_y2);
        vec_res3 = _mm256_fmadd_pd(vec_a3,vec_x3,vec_y3);
        vec_res4 = _mm256_fmadd_pd(vec_a4,vec_x4,vec_y4);
        vec_res5 = _mm256_fmadd_pd(vec_a5,vec_x5,vec_y5);
        vec_res6 = _mm256_fmadd_pd(vec_a6,vec_x6,vec_y6);
        vec_res7 = _mm256_fmadd_pd(vec_a7,vec_x7,vec_y7);
        vec_res8 = _mm256_fmadd_pd(vec_a8,vec_x8,vec_y8);


        _mm256_storeu_pd(y+i   ,vec_res1);
        _mm256_storeu_pd(y+i+4 ,vec_res2);
        _mm256_storeu_pd(y+i+8 ,vec_res3);
        _mm256_storeu_pd(y+i+12,vec_res4);
        _mm256_storeu_pd(y+i+16,vec_res5);
        _mm256_storeu_pd(y+i+20,vec_res6);
        _mm256_storeu_pd(y+i+24,vec_res7);
        _mm256_storeu_pd(y+i+28,vec_res8);
    }
    for(; i < n; ++i ) {
        y [i] += x[i] * a;
    }
}
//z <- a * x + y
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
    auto peel = (unsigned long) x & 0x1f;
    if (peel != 0){
        peel = ( 32-peel) *d_size_inv;
        for(; i < peel ;++i){
            z[i] = x[i] * a + y[i];
        }
    }
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


// y <- a * x + y; returns max |y[i]|
/*
double axpymax(const unsigned int n, const double a, const double *x, double *y) {
	double max_abs_val = 0;
	for(unsigned i = 0 ; i < n; ++i ) {
		y [i] += x[i] * a;
		const double abs_val = std::abs(y[i]);
		if (max_abs_val < abs_val) max_abs_val = abs_val;
	}
	return max_abs_val;
}
*/
double axpymax(const unsigned int n, const double a, const double *x, double *y) {
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
    const __m256d  vec_a1 = _mm256_set1_pd(a);
    const __m256d  vec_a2 = _mm256_set1_pd(a);
    const __m256d  vec_a3 = _mm256_set1_pd(a);
    const __m256d  vec_a4 = _mm256_set1_pd(a);
    const __m256d  vec_a5 = _mm256_set1_pd(a);
    const __m256d  vec_a6 = _mm256_set1_pd(a);
    const __m256d  vec_a7 = _mm256_set1_pd(a);
    const __m256d  vec_a8 = _mm256_set1_pd(a);
    __m256d vec_max1 = _mm256_setzero_pd();
    __m256d vec_max2 = _mm256_setzero_pd();
    __m256d vec_max3 = _mm256_setzero_pd();
    __m256d vec_max4 = _mm256_setzero_pd();
    __m256d vec_max5 = _mm256_setzero_pd();
    __m256d vec_max6 = _mm256_setzero_pd();
    __m256d vec_max7 = _mm256_setzero_pd();
    __m256d vec_max8 = _mm256_setzero_pd();

    // peel loop for alligned loads for a ! b should be handled too.
    auto peel = (unsigned long) x & 0x1f;
    if (peel != 0){
        peel = ( 32-peel) *d_size_inv;
        for(; i < peel ;++i){
            y [i] += x[i] * a;
            const double square_val = y[i]* y[i];
            if (max_abs_val < square_val) max_abs_val = square_val;
        }
    }
    for(; i < n-32; i +=32 ) {
        vec_x1      = _mm256_load_pd(x+i);
        vec_x2      = _mm256_load_pd(x+i+4);
        vec_x3      = _mm256_load_pd(x+i+8);
        vec_x4      = _mm256_load_pd(x+i+12);
        vec_x5      = _mm256_load_pd(x+i+16);
        vec_x6      = _mm256_load_pd(x+i+20);
        vec_x7      = _mm256_load_pd(x+i+24);
        vec_x8      = _mm256_load_pd(x+i+28);
        vec_y1      = _mm256_load_pd(y+i);
        vec_y2      = _mm256_load_pd(y+i+4);
        vec_y3      = _mm256_load_pd(y+i+8);
        vec_y4      = _mm256_load_pd(y+i+12);
        vec_y5      = _mm256_load_pd(y+i+16);
        vec_y6      = _mm256_load_pd(y+i+20);
        vec_y7      = _mm256_load_pd(y+i+24);
        vec_y8      = _mm256_load_pd(y+i+28);

        // not sure if we get aliasing issues here
        vec_res1 = _mm256_fmadd_pd(vec_a1,vec_x1,vec_y1);
        vec_res2 = _mm256_fmadd_pd(vec_a2,vec_x2,vec_y2);
        vec_res3 = _mm256_fmadd_pd(vec_a3,vec_x3,vec_y3);
        vec_res4 = _mm256_fmadd_pd(vec_a4,vec_x4,vec_y4);
        vec_res5 = _mm256_fmadd_pd(vec_a5,vec_x5,vec_y5);
        vec_res6 = _mm256_fmadd_pd(vec_a6,vec_x6,vec_y6);
        vec_res7 = _mm256_fmadd_pd(vec_a7,vec_x7,vec_y7);
        vec_res8 = _mm256_fmadd_pd(vec_a8,vec_x8,vec_y8);


        _mm256_storeu_pd(y+i   ,vec_res1);
        _mm256_storeu_pd(y+i+4 ,vec_res2);
        _mm256_storeu_pd(y+i+8 ,vec_res3);
        _mm256_storeu_pd(y+i+12,vec_res4);
        _mm256_storeu_pd(y+i+16,vec_res5);
        _mm256_storeu_pd(y+i+20,vec_res6);
        _mm256_storeu_pd(y+i+24,vec_res7);
        _mm256_storeu_pd(y+i+28,vec_res8);

        vec_res1 = _mm256_mul_pd(vec_res1,vec_res1);
        vec_res2 = _mm256_mul_pd(vec_res2,vec_res2);
        vec_res3 = _mm256_mul_pd(vec_res3,vec_res3);
        vec_res4 = _mm256_mul_pd(vec_res4,vec_res4);
        vec_res5 = _mm256_mul_pd(vec_res5,vec_res5);
        vec_res6 = _mm256_mul_pd(vec_res6,vec_res6);
        vec_res7 = _mm256_mul_pd(vec_res7,vec_res7);
        vec_res8 = _mm256_mul_pd(vec_res8,vec_res8);


        vec_max1 = _mm256_max_pd(vec_res1,vec_max1);
        vec_max2 = _mm256_max_pd(vec_res2,vec_max2);
        vec_max3 = _mm256_max_pd(vec_res3,vec_max3);
        vec_max4 = _mm256_max_pd(vec_res4,vec_max4);
        vec_max5 = _mm256_max_pd(vec_res5,vec_max5);
        vec_max6 = _mm256_max_pd(vec_res6,vec_max6);
        vec_max7 = _mm256_max_pd(vec_res7,vec_max7);
        vec_max8 = _mm256_max_pd(vec_res8,vec_max8);

    }
   for(; i < n; ++i ) {
       y [i] += x[i] * a;
       const double square_val = y[i]*y[i];
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
   vec_max2 = _mm256_permute_pd(vec_max1,0xb1);
   vec_max1 = _mm256_max_pd(vec_max2,vec_max1);
   // sanity check: here element 0 and 1 aswell as 2 and 3 should be the same

    // maybe add allignement
    double sol0 [4];
    _mm256_storeu_pd(sol0,vec_max1);

    //potential optimisation :
    //const __m256i mask = {1,0,1,0}; // maybe 256
    //_mm256_maskstore_pd(sol0,0xff00ff00,vec_max1);

    return std::sqrt(std::max(std::max(sol0[0],sol0[2]),max_abs_val));

}

// z <- a * x + y; returns max |z[i]|
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
    auto peel = (unsigned long) x & 0x1f;
    if (peel != 0) {
        peel = (32 - peel) * d_size_inv;
        for (; i < peel; ++i) {
            z[i] = y[i] + x[i] * a;
            const double square_val = z[i] * z[i];
            if (max_abs_val < square_val) max_abs_val = square_val;
        }
    }
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


// returns max |x[i]| for i in 0:n-1
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
    unsigned peel = (unsigned long) x & 0x1f;
    if (peel != 0) {
        peel = (32 - peel) * d_size_inv;
        for (; i < peel; ++i) {
            const double square_val = x[i]*x[i];
            if (max_abs_val < square_val) max_abs_val = square_val;
        }
    }
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

ICConjugateGradientSolver::ICConjugateGradientSolver(unsigned max_steps, const Mac3d& grid)
   :
       grid{grid},
       n_cells_x{grid.get_num_cells_x()}, n_cells_y{grid.get_num_cells_y()}, n_cells_z{grid.get_num_cells_z()},
       tau{0.97},
       num_cells{n_cells_x * n_cells_y * n_cells_z},
       max_steps(max_steps)
{
   step = 0;
   unsigned chunk_size = 32;
   q = new (std::align_val_t(chunk_size)) double [num_cells];
   r = new (std::align_val_t(chunk_size)) double [num_cells];
   z = new (std::align_val_t(chunk_size)) double [num_cells];
   s = new (std::align_val_t(chunk_size)) double [num_cells];
   precon_diag = new (std::align_val_t(chunk_size)) double [num_cells];
   A_diag = new (std::align_val_t(chunk_size)) double [num_cells];

   // load A_diag
   unsigned cellidx = 0;
   for (unsigned k = 0; k < n_cells_z; k++) {
       for (unsigned j = 0; j < n_cells_y; j++) {
           for (unsigned i = 0; i < n_cells_x; i++, cellidx++) {
               A_diag[cellidx] = grid.A_diag_val[i + j*stride_y + k*stride_z];
           }
       }
   }
}

// apply the preconditioner (L L^T)^-1 by solving Lq = d and Lp = q
void ICConjugateGradientSolver::applyPreconditioner(const double *r, double *z) const {
   for (unsigned k = 0; k < n_cells_z; k++) {
       for (unsigned j = 0; j < n_cells_y; j++) {
           for (unsigned i = 0; i < n_cells_x; i++) {
               const unsigned cellidx = i + j*stride_y + k*stride_z;
               if (not grid.pfluid_[cellidx]) {
                   q[cellidx] = 0;
                   continue;
               }
               double t = r[cellidx];
               if (i > 0 && grid.pfluid_[cellidx - stride_x]) t += precon_diag[cellidx - stride_x] * q[cellidx - stride_x];
               if (j > 0 && grid.pfluid_[cellidx - stride_y]) t += precon_diag[cellidx - stride_y] * q[cellidx - stride_y];
               if (k > 0 && grid.pfluid_[cellidx - stride_z]) t += precon_diag[cellidx - stride_z] * q[cellidx - stride_z];
               q[cellidx] = t * precon_diag[cellidx];
           }
       }
   }
   for (int k = n_cells_z-1; k >= 0; k--) {
       for (int j = n_cells_y-1; j >= 0; j--) {
           for (int i = n_cells_x-1; i >= 0; i--) {
               const unsigned cellidx = i + j*stride_y + k*stride_z;
               if (not grid.pfluid_[cellidx]) {
                   z[cellidx] = 0;
                   continue;
               }
               double t = q[cellidx];
               if (i + 1 < (int) n_cells_x && grid.pfluid_[cellidx + stride_x]) t += precon_diag[cellidx] * z[cellidx + stride_x];
               if (j + 1 < (int) n_cells_y && grid.pfluid_[cellidx + stride_y]) t += precon_diag[cellidx] * z[cellidx + stride_y];
               if (k + 1 < (int) n_cells_z && grid.pfluid_[cellidx + stride_z]) t += precon_diag[cellidx] * z[cellidx + stride_z];
               z[cellidx] = t * precon_diag[cellidx];
           }
       }
   }
}

void ICConjugateGradientSolver::computePreconDiag() {
   for (unsigned k = 0; k < n_cells_z; k++) {
       for (unsigned j = 0; j < n_cells_y; j++) {
           for (unsigned i = 0; i < n_cells_x; i++) {
               const unsigned cellidx = i + j*stride_y + k*stride_z;
               if (not grid.pfluid_[cellidx]) {
                   precon_diag[cellidx] = 0;
                   continue;
               }
               double e = A_diag[cellidx];
               if (i > 0 && grid.pfluid_[cellidx - stride_x]) e -= (precon_diag[cellidx - stride_x] * precon_diag[cellidx - stride_x]);
               if (j > 0 && grid.pfluid_[cellidx - stride_y]) e -= (precon_diag[cellidx - stride_y] * precon_diag[cellidx - stride_y]);
               if (k > 0 && grid.pfluid_[cellidx - stride_z]) e -= (precon_diag[cellidx - stride_z] * precon_diag[cellidx - stride_z]);
               precon_diag[cellidx] = 1 / std::sqrt(e + 1e-30);
           }
       }
   }
}

// apply the matrix A: y <- A b
// element order: k-j-i
void ICConjugateGradientSolver::applyA(const double *b, double *y) const{
   for (unsigned k = 0; k < n_cells_z; k++) {
       for (unsigned j = 0; j < n_cells_y; j++) {
           for (unsigned i = 0; i < n_cells_x; i++) {
               const unsigned cellidx = i + j*stride_y + k*stride_z;
               double diag_val = grid.A_diag_val[cellidx];
               const double b_cell = b[cellidx];

               if (i == 0 && j == 0 && k == 0) y[cellidx] = 0;
               y[cellidx] += diag_val * b_cell;

               // avoid filling y with zeroes beforehand
               if (j == 0 && k == 0 && i+1 < n_cells_x)  y[cellidx+stride_x] = 0;
               if (k == 0 && j+1 < n_cells_y)  y[cellidx+stride_y] = 0;
               if (k+1 < n_cells_z) y[cellidx+stride_z] = 0;
               // Compute off-diagonal entries
               if (grid.pfluid_[cellidx]){
                   // x-adjacent cells
                   if (i+1 < n_cells_x && grid.pfluid_[cellidx+stride_x]){
                       y[cellidx] 			-= b[cellidx+stride_x];
                       y[cellidx+stride_x] -= b_cell;
                   }

                   // y-adjacent cells
                   if (j+1 < n_cells_y && grid.pfluid_[cellidx+stride_y]){
                       y[cellidx] 			-= b[cellidx+stride_y];
                       y[cellidx+stride_y] -= b_cell;
                   }

                   // z-adjacent cells
                   if (k+1 < n_cells_z && grid.pfluid_[cellidx+stride_z]){
                       y[cellidx] 			-= b[cellidx+stride_z];
                       y[cellidx+stride_z] -= b_cell;
                   }
               }
           }
       }
   }
}

void checknan(const double* array, int len, std::string array_name="array") {
   int count_left = 20;
   for (int i = 0; i<len; i++) {
       if(std::isnan(array[i])) {
           std::cout
               << "WARNING: NaN found in " << array_name <<
               " at position " << i << "!" << std::endl;
           if (--count_left == 0)
               std::cout << "[further occurrences ignored]" << std::endl;
       }
   }
}

void print_array_head(const double* array, std::string prefix="", unsigned number=20) {
   std::cout << prefix;
   for (unsigned i=0; i < number; i++) std::cout << array[i] << ' ';
   std::cout << std::endl;
}


void ICConjugateGradientSolver::solve(const double* rhs, double* p) {
   // initialize initial guess and residual
   // catch zero rhs early
   double max_residual_modulus = xmax(num_cells, rhs);
   if (max_residual_modulus < thresh) {
       std::fill(p,p+num_cells,0);
       return;
   }

   computePreconDiag();
   // s = M⁻¹ r
   applyPreconditioner(rhs, s);

   // ρ = <r,s>
   double rho = dot(rhs,s,num_cells);

   for(unsigned step = 0; step < max_steps; step++){
       applyA(s, z);
       const double dots = dot(z,s,num_cells);
       const double alpha = rho / dots;

       double max_abs_val;
       if (step == 0) {
           // on the first step initialize p
           // p <- α s
           axy(num_cells, alpha, s, p);

           // r <- (-α z + rhs)
           max_abs_val = axpyzmax(num_cells, -alpha, z, rhs, r);
       }
       else {
           // p <- α s
           axpy(num_cells, alpha, s, p);

           // r <- (-α z + r)
           max_abs_val = axpymax(num_cells, -alpha, z, r);
       }

       // std::cout << "Max abs val: " << max_abs_val << std::endl;
       if (max_abs_val < thresh) {
           // std::cout << "Number of steps: " << step + 1 << std::endl;
           return;
       }
       else if (step + 1 == max_steps) {
           std::cout << "WARNING: Failed to find a solution in " << step + 1 << " steps (|r|=" << max_abs_val <<")!" << std::endl;
       }

       // z = M⁻¹ r
       applyPreconditioner(r, z);
       const double rho_new = dot(z, r, num_cells);
       const double beta = rho_new / rho;
       rho = rho_new;
       //Bug potential: aliasing
       axpyz(num_cells, beta, s, z, s);
   }
}

SparseMat::SparseMat(unsigned a, unsigned b): v(a), r(b){
   values = new (std::align_val_t(32)) double [v];
   col_idx = new (std::align_val_t(32)) unsigned [v];
   row_idx = new (std::align_val_t(32)) unsigned [r];
}
SparseMat::~SparseMat(){
   delete [] values;
   delete [] col_idx;
   delete [] row_idx;
}


ICConjugateGradientSolver::~ICConjugateGradientSolver() {
   delete [] q;
   delete [] z;
   delete [] s;
   delete [] precon_diag;
   delete [] A_diag;

}


