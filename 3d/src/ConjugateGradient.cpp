#include "ConjugateGradient.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

#include "kernels/dot.h"
#include "kernels/axy.h"
#include "kernels/axpy.h"
#include "kernels/axpyz.h"
#include "kernels/axpymax.h"
#include "kernels/axpyzmax.h"
#include "kernels/xmax.h"

#define square(X)	((X) * (X))

// ********* Kernels **********

// variables for optimisation
// const unsigned d_size_inv = 1/sizeof(double);



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


