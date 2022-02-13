# Benchmarks

All benchmarks have fixed random seed for deterministic execution.
All benchmarks have "show wireframe" enabled to force computation of mesh at each step, but "Export mesh" disabled.

Note: I'm reporting indicative RAM usage and runtime for each benchmark. These measurements were taken before extensive optimisations at tag `v1.1`. Of course these numbers are purely indicative and subject to change as we optimize and depend on the machine running the benchmarks.

- **Case 1:** cube of fluid at the center of the domain.
  - `benchmark-1-0`: 10x10x10 grid; 512 particles; 1000 time steps (~1.6MB RAM, <1min). 
  - `benchmark-1-1`: 20x20x20 grid; 8,000 particles; 1000 time steps (~4MB RAM, ~2mins).
  - `benchmark-1-2`: 40x40x40 grid; 64,000 particles; 1000 time steps (~26MB RAM, ~30mins).
  - `benchmark-1-3`: 80x80x80 grid; 512,000 particles; 500 time steps (~200MB RAM, ~2hrs).
  - `benchmark-1-4`: 160x160x160 grid; 4,096,000 particles; 40 time steps (~1.1-1.4GB RAM, ~1hr).
- **Case 2:** dam break.
  - `benchmark-2-0`: 30x10x5 grid; 1,920 particles; 1500 time steps (~1.2MB RAM, <1min).
  - `benchmark-2-1`: 60x20x10 grid; 1,4080 particles; 1500 time steps (~5.5MB RAM, ~5mins).
  - `benchmark-2-2`: 120x40x20 grid; 10,4160 particles; 750 time steps (~40MB RAM, ~40mins).
  - `benchmark-2-3`: 240x80x40 grid; 800,320 particles; 750 time steps (~300MB RAM, ~3.5hrs).
  - `benchmark-2-4`: 480x160x80 grid; 6,272,640 particles; 20 time steps (~1.4-2.1GB  RAM, ~1hr).

## Running benchmarks

Benchmarks are subdivided in two categories: "fast" benchmarks, which should take less than a minute each,
and "slow" benchmarks, which may take several hours each.

To try out just the "fast" benchmarks, run `make fast`.
If you want to run the full suite, run `make all`, but beware that it may take hours.

The default build directory is `../build/`. It is assumed that the executable `watersim-cli` is present there.
To run with a custom build directory, run `make build-directory=../foo/bar/ [all|fast|slow]`.

You can run a single benchmark by specifying its config file name, e.g. `make benchmark-1-0.json`.

To remove all benchmark data, use `make clean`.

## Benchmarking checkpoints

Below is a list of git tags which represent "benchmarking checkpoints" in development.
These were used to perform full benchmarks in order to measure performance improvements. See ASL-FLIP report for details.

- `v1.1`: baseline benchmark. No optimizations performed, only restructuring and instrumentation.
- `v1.2`: Optimizations to level-set computation.
- `v1.3`: New and improved Particles data structure.
- `v1.4`: Complete rewrite of interpolation methods.
- `v1.5`: Optimized particle-to-grid, forces and BCs.
- `v1.6`: Custom Conjugate Gradient solver.
