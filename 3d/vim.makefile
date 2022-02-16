
# This file allows for convenient compilation and execution
# of the project from inside vim (or any other IDE for that matter)

build-directory ?= .vim-build
jobs=8

all: _vim_run

# Commands for Vim to compile and run the project
.PHONY: _vim_run _vim_build _vim_clean _vim_test
_vim_run: _vim_build
	cd $(build-directory)
	./watersim-gui

_vim_build: watersim-gui

_vim_test: watersim-tests

_vim_clean: clean
	rm -rf $(build-directory)

$(build-directory)/Makefile:
	mkdir -p $(build-directory)
	cd $(build-directory)
	cmake ..

binaries = watersim-cli watersim-gui watersim-tests
$(binaries): $(build-directory)/Makefile .FORCE
	cd $(build-directory)
	make -j$(jobs) $@

.PHONY: .ONESHELL .FORCE
.ONESHELL:


