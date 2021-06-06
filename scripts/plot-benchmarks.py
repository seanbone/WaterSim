#!/usr/bin/env python3

"""
This script plots benchmark data for WaterSim. Run with -h for usage instructions.

The script scans the search directory for benchmark folders. Each benchmark folder should contain the output of
running the full benchmark suite at a given optimization point (git tag).
Recommended directory structure:

|-- benchmarks -> directory to be scanned
    |-- benchmarks-full-v1.1 -> output of benchmarks at tag v1.1
    |   |-- benchmark-1-0
    |   |-- benchmark-1-1
    |   |-- [...]
    |   `-- benchmark-2-4
    |-- benchmarks-full-v1.2
    |   |-- benchmark-1-0
    |   |-- benchmark-1-1
    |   |-- [...]
    |   `-- benchmark-2-4
    `-- [...]
"""
import os
import argparse
import matplotlib.pyplot as plt
from cycler import cycler
import numpy as np


def get_directories(path):
    """ Generates a list of directories at given path. """
    filenames = os.listdir(os.path.abspath(path))

    dirs = []
    for filename in filenames:
        if os.path.isdir(os.path.join(os.path.abspath(path), filename)):
            dirs.append(filename)

    dirs.sort()
    return dirs


def get_benchmark_names(num_cases=2, num_per_case=5):
    """ Generates the list of benchmark names. """
    names = []
    for i in range(1, num_cases+1):
        for j in range(num_per_case):
            names.append(f"benchmark-{i}-{j}")

    return names


def collect_all_timings(path, prefix, labels, substeps):
    """
    Collect all timing data.
    Example of access to return dictionary: data['v1.1']['benchmark-1-0']['FLIP']['mean']
    """
    data = {}
    for label in labels:
        data[label] = {}
        batch_path = os.path.join(path, f"{prefix}{label}")
        benchmark_names = get_benchmark_names()

        for b in benchmark_names:
            filename = os.path.join(batch_path, b, "timing_info.txt")
            print(f"Reading {filename}...")

            filename = os.path.abspath(filename)
            if not os.path.isfile(filename):
                exit('Error: could not open file.')

            with open(filename) as f:
                lines = f.read().splitlines()

            data[label][b] = {}
            for line in lines:
                for substep in substeps:
                    if line.find(substep) > -1:
                        data[label][b][substep] = {}
                        vals = line.split()[1].split('/')
                        data[label][b][substep]["mean"] = float(vals[0])
                        data[label][b][substep]["median"] = float(vals[1])
                        data[label][b][substep]["stdev"] = float(vals[2])
        
    return data


def collect_all_flops(path, prefix, labels, substeps):
    """
    Collect all flops data.
    Example of access to return dictionary: data['v1.1']['benchmark-1-0']['FLIP']['mean']
    """
    data = {}
    
    for label in labels:
        data[label] = {}
        batch_path = os.path.join(path, f"{prefix}{label}")
        benchmark_names = get_benchmark_names()

        for b in benchmark_names:
            filename = os.path.join(batch_path, b, "cost_analysis.txt")
            print(f"Reading {filename}...")

            filename = os.path.abspath(filename)
            if not os.path.isfile(filename):
                exit('Error: could not open file.')

            with open(filename) as f:
                lines = f.read().splitlines()

            data[label][b] = {}
            for line in lines:
                for substep in substeps:
                    if line.find(substep) > -1:
                        data[label][b][substep] = {}
                        vals = line.split()
                        data[label][b][substep]['adds'] = float(vals[1])
                        data[label][b][substep]['muls'] = float(vals[2])
                        data[label][b][substep]['divs'] = float(vals[3])
                        data[label][b][substep]['read'] = float(vals[4])
        
    return data


def plot_histogram(benchmark, tags, substeps, all_data, output, show=False, title=None):
    """ Generate a histogram of runtimes for a specific benchmark across all optimization tags. """
    if title is None:
        title = f"Runtimes for benchmark '{benchmark}'"

    my_data = {}
    for substep in substeps:
        tmp = []
        for tag in tags:
            tmp.append(all_data[tag][benchmark][substep]['mean'])
        my_data[substep] = np.array(tmp)

    num_tags = len(tags)
    bottoms = np.zeros(num_tags)
    indices = range(num_tags)

    fig, ax = plt.subplots()

    for substep in substeps:
        ax.bar(indices, my_data[substep], bottom=bottoms, label=substep, tick_label=tags)
        bottoms += my_data[substep]

    ax.legend()
    ax.set_ylabel("Average runtime [cycles]")
    ax.set_xlabel("Optimization stage")
    ax.set_title(title)
    if show:
        plt.show()
    else:
        filename = output.replace("%", benchmark)
        print(f"Writing '{filename}'...")
        fig.savefig(filename, dpi=300)
    plt.close(fig)


def plot_perf(tag, substeps, all_data, all_flops, output, show=False, title=None):
    if title is None:
        title = f"Runtime plot for optimization stage '{tag}'"

    benchmarks = get_benchmark_names(1)
    problem_dimensions = [10**3, 20**3, 40**3, 80**3, 160**3]

    my_data = {}
    for substep in substeps:
        tmp = []
        for benchmark in benchmarks:
            if substep != "compute_mesh": tmp.append((all_flops[tag][benchmark][substep]['adds'] + all_flops[tag][benchmark][substep]['muls']) / all_data[tag][benchmark][substep]['mean'])
            else:                         tmp.append(0.)
            # tmp.append(all_data[tag][benchmark][substep]['mean'])
        my_data[substep] = np.array(tmp)

    fig, ax = plt.subplots()

    my_cycler = (cycler(marker="ovsXD") + (cycler(color=['tab:red', 'tab:cyan', 'tab:blue', 'tab:orange', 'tab:green'])+
                 cycler(linestyle=['-', '--', ':', '-.', '--'])))
    ax.set_prop_cycle(my_cycler)


    ax.set_ylabel("Average performance [cycles]")
    ax.set_xlabel("Number of cells [-]")
    for substep in substeps:
        ax.plot(problem_dimensions, my_data[substep], label=substep)

    ax.legend()
    ax.set_yscale("log", basey=10)
    ax.set_xscale("log", basex=2)
    ax.set_title(title)
    ax.set_facecolor('#cccccc')
    ax.grid(axis='both', which='major', zorder=-1, color='white')
    if show:
        plt.show()
    else:
        filename = output.replace("%", f"perf-{tag}")
        print(f"Writing '{filename}'...")
        fig.savefig(filename, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    path = "benchmarks/"
    prefix = "benchmarks-full-"
    substeps = ["particle_to_grid", "apply_pressure_correction",
                "grid_to_particle", "advance_particles", "compute_mesh"]
    output = "%.pdf"

    parser = argparse.ArgumentParser(description="Generate plots of the WaterSim benchmarks.")

    parser.add_argument("-a", "--all", action="store_true", help="Generate all possible plots.")
    parser.add_argument("-b", "--benchmark", nargs=1, action="append", help="Generate a histogram comparing a specific benchmark across each tag. Can be called multiple times for more than one histogram.")
    parser.add_argument("-d", "--dir", nargs=1, default=[path], help=f"Directory to search for benchmarks. Default %(default)s.")
    parser.add_argument("-o", "--output", nargs=1, default=[output], help=f"Output file name pattern. The character '%%' will be replaced with the file name. Default: %(default)s")
    parser.add_argument("--prefix", nargs=1, default=[prefix], help=f"Prefix of benchmark tag directories. Default %(default)s.")
    parser.add_argument("-s", "--substep", nargs=1, action="append", help="Specify a FLIP substep to include in plots. Can be called multiple times for more than one substep.")
    parser.add_argument("--show", action="store_true", help="Show the plots instead of writing to file(s).")
    parser.add_argument("-t", "--tag", nargs=1, action="append", help="Generate a performance plot (problem size vs runtime) for the specified tag. Can be called multiple times for more than one plot.")

    args = parser.parse_args()
    #print(args)
    path = args.dir[0]
    prefix = args.prefix[0]
    substeps = [s[0] for s in args.substep] if args.substep is not None else substeps
    output = args.output[0]
    show = args.show

    if not args.all and args.benchmark is None and args.tag is None:
        exit("Warning: no plots to be generated!\nRun with '-h' for help.")

    dirs = get_directories(path)

    all_tags = []
    for d in dirs:
        if d.find(prefix) == 0:
            all_tags.append(d.replace(prefix, "", 1))

    print(all_tags, args.tag)

    all_data  = collect_all_timings(path, prefix, all_tags, substeps)
    all_flops = collect_all_flops(path, prefix, all_tags, substeps)

    benchmarks_hist_plots = None
    if args.all:
        benchmarks_hist_plots = get_benchmark_names()
    elif args.benchmark is not None:
        benchmarks_hist_plots = [b[0] for b in args.benchmark]

    if benchmarks_hist_plots is not None:
        for b in benchmarks_hist_plots:
            plot_histogram(b, all_tags, substeps, all_data, output, show)

    benchmark_perf_plots = None
    if args.all:
        benchmark_perf_plots = all_tags
    elif args.tag is not None:
        benchmark_perf_plots = [t[0] for t in args.tag]

    if benchmark_perf_plots is not None:
        for tag in benchmark_perf_plots:
            plot_perf(tag, substeps, all_data, all_flops, output, show)
