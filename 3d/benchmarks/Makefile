# A Makefile for easy automagical execution of benchmarks.
# If you don't have much time, just run `make fast`
# If you want to run the full suite, run `make all`, but beware that it may take hours.
# To run with custom build directory, run `make build-directory=..foo/bar/`

# Note: executable watersim-cli is assumed to exist in build directory
build-directory ?= ../build/

post-script ?= ../../scripts/load_timing_info.py

# List of all benchmarks to run. Divided into "fast" and "slow" benchmarks
fast-benchmarks = benchmark-1-0 benchmark-2-0
slow-benchmarks = benchmark-1-1 benchmark-1-2 benchmark-1-3 benchmark-1-4 benchmark-2-1 benchmark-2-2 benchmark-2-3 benchmark-2-4

fast-benchmark-files = $(addsuffix .json, $(fast-benchmarks))
slow-benchmark-files = $(addsuffix .json, $(slow-benchmarks))

.PHONY: warning fast slow all
warning:
	echo "Note: only running 'fast' benchmarks. To run all benchmarks, use 'make all'."
	make fast

fast: $(fast-benchmark-files)

slow: $(slow-benchmark-files)

all: fast slow

$(fast-benchmark-files) $(slow-benchmark-files): .FORCE
	echo -n "Running benchmark '$(basename $@)'... "
	mkdir -p $(basename $@)
	cd $(basename $@)
	../$(build-directory)watersim-cli -y ../$@ > $(basename $@).log
	python3 ../$(post-script) timings.json > timing_info.txt
	echo "Done."

.PHONY: .FORCE
.FORCE:

# .ONESHELL and .SILENT are special targets, see
# https://www.gnu.org/software/make/manual/html_node/Special-Targets.html
.ONESHELL:

ifndef VERBOSE
.SILENT:
endif

.PHONY: clean
clean:
	rm -rf $(fast-benchmarks) $(slow-benchmarks)