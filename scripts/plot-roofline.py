import os
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np


# Peak performance [flops/cycle]
P     = 4.  # 2 FMAs
P_vec = 16. # 256-bit registers

# Peak bandwidth [GB/s]
B_theoretical = 19.2    # DDR clock rate * #channels * bits per channel per cycle * 0.125 B/bit
B_stream      = 10.     # STREAM benchmark (https://github.com/jeffhammond/STREAM)

# Processor frequency [GHz]
freq = 1.8

# Peak bandwidth [B/cycle]
Bc_theoretical = B_theoretical/freq
Bc_stream      = B_stream/freq

# Ridge point [flops/B]
ridge_point_theoretical     = P/Bc_theoretical
ridge_point_stream          = P/Bc_stream
ridge_point_theoretical_vec = P_vec/Bc_theoretical
ridge_point_stream_vec      = P_vec/Bc_stream

# Fonts
font            = {"family": "DejaVu Serif"}
annotation_font = {"family": "DejaVu Serif", "size": 8}
code_font       = {"family": "Courier New"}
legend_settings = {"family": "DejaVu Serif", "size": 8}

# Read flops and memops values
n_sections   = 4
n_operations = 4
ops    = np.empty([n_sections, n_operations])
cycles = np.empty([n_sections])
for b in range(2):      # Benchmark type
    for s in range(5):  # Benchmark size
        
        # TODO: Expand so that here we have a loop over the tags

        ops_filename    = "../3d/benchmarks/benchmark-" + str(b+1) + "-" + str(s) + "/cost_analysis.txt"
        cycles_filename = "../3d/benchmarks/benchmark-" + str(b+1) + "-" + str(s) + "/timing_info.txt"
        
        if os.path.isfile(ops_filename) and os.path.isfile(cycles_filename):
            
            ops_file    = open(ops_filename,    "r")
            cycles_file = open(cycles_filename, "r")

            counter = 0
            for line in ops_file:   # Section type
                counter += 1
                if counter > 4:
                    tmp = line.split()
                    
                    for i in range(1, 5):   # Operation type
                        ops[counter-5, i-1] = float(tmp[i])

            counter = 0

            for line in cycles_file:    # Section type
                counter += 1
                if counter == 5:
                    cycles[counter-5] = float(line.split()[1].split("/")[0])
                elif counter >= 8 and counter <= 10:
                    cycles[counter-7] = float(line.split()[1].split("/")[0])
            
            ops_file.close()
            cycles_file.close()

            ptg_perf = (ops[0, 0] + ops[0, 1]) / cycles[0]
            prc_perf = (ops[1, 0] + ops[1, 1]) / cycles[1]
            gtp_perf = (ops[2, 0] + ops[2, 1]) / cycles[2]
            adp_perf = (ops[3, 0] + ops[3, 1]) / cycles[3]

            ptg_oint = (ops[0, 0] + ops[0, 1]) / ops[0, 3]
            prc_oint = (ops[1, 0] + ops[1, 1]) / ops[1, 3]
            gtp_oint = (ops[2, 0] + ops[2, 1]) / ops[2, 3]
            adp_oint = (ops[3, 0] + ops[3, 1]) / ops[3, 3]

            # Plot the data
            fig, ax = plt.subplots()

            # Theoretical rooflines
            ax.plot([2.**(-5.),                   ridge_point_theoretical    ], [(2**(-5)) * Bc_theoretical, P    ], color="#177245", linestyle="-" ) # Memory Bound Scalar
            ax.plot([ridge_point_theoretical,     ridge_point_theoretical_vec], [P,                          P_vec], color="#03c03c", linestyle="--") # Memory Bound Vectorized
            ax.plot([ridge_point_theoretical,     2.**5.                     ], [P,                          P    ], color="#177245", linestyle="-" ) # Compute Bound Scalar
            ax.plot([ridge_point_theoretical_vec, 2.**5.                     ], [P_vec,                      P_vec], color="#03c03c", linestyle="--") # Compute Bound Vectorized

            # STREAM rooflines
            ax.plot([2.**(-5.),              ridge_point_stream    ], [(2**(-5)) * Bc_stream, P    ], color="#b04d4f", linestyle="-" ) # Memory Bound Scalar
            ax.plot([ridge_point_stream,     ridge_point_stream_vec], [P,                     P_vec], color="#e8792b", linestyle="--") # Memory Bound Vectorized
            ax.plot([ridge_point_stream,     2.**5.                ], [P,                     P    ], color="#b04d4f", linestyle="-" ) # Compute Bound Scalar
            ax.plot([ridge_point_stream_vec, 2.**5.                ], [P_vec,                 P_vec], color="#e8792b", linestyle="--") # Compute Bound Vectorized

            # Annotations
            ax.annotate("Scalar roofline",       (1.333 * ridge_point_theoretical_vec, 1.0625 * P    ),              color="darkslategrey", **annotation_font)
            ax.annotate("Vectorized roofline",   (1.333 * ridge_point_theoretical_vec, 1.0625 * P_vec),              color="darkslategrey", **annotation_font)
            ax.annotate("STREAM bandwidth",      (0.25  * ridge_point_stream,          0.3    * P    ), rotation=38, color="darkslategrey", **annotation_font)
            ax.annotate("Theoretical bandwidth", (0.25  * ridge_point_theoretical,     0.3    * P    ), rotation=38, color="darkslategrey", **annotation_font)
            
            # Plot performances
            ax.scatter(ptg_oint, ptg_perf, color="#148629", marker="x", zorder=3, label="Particle-to-grid"   )
            ax.scatter(prc_oint, prc_perf, color="#226095", marker="v", zorder=3, label="Pressure correction")
            ax.scatter(gtp_oint, gtp_perf, color="#d96465", marker="o", zorder=3, label="Grid-to-particle"   )
            ax.scatter(adp_oint, adp_perf, color="#8552ad", marker="+", zorder=3, label="Advance particles"  )

            # Set legend
            legend = ax.legend(loc=2, ncol=2, prop=annotation_font)
            frame = legend.get_frame()
            frame.set_facecolor("none")
            frame.set_edgecolor("darkslategrey")
            frame.set_boxstyle("square", pad=0.)
            
            # Set background color
            ax.set_facecolor("lightgrey")

            # Set border of the image to not be visible
            ax.spines["left"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines['top'].set_visible(False)

            # Set aspect of the plot: auto --> fill the position rectangle with data
            ax.set_aspect("auto")

            # Set the axes scale
            ax.set_xscale("log", base=2)
            ax.set_yscale("log", base=2)

            # Show full numbers on axes
            #ax.xaxis.set_major_formatter(ScalarFormatter())
            #ax.yaxis.set_major_formatter(ScalarFormatter())

            # Set title and labels of the axes
            ax.set_title("Roofline plot", loc="left", pad=14., fontweight="bold", **font)
            ax.set_xlabel("Operational intensity [flops/byte]", **font)
            ax.set_ylabel("Performance [flops/cycle]", rotation="horizontal", loc="bottom", **font)

            # Zoom
            ax.set_xlim(2.**(-4.), 2.**3.)
            ax.set_ylim(2.**(-2.), 2.**5.)

            ax.yaxis.grid(True, which='major', linestyle='-', color="white")
            ax.yaxis.set_tick_params(length=0)

            ax.yaxis.set_label_coords(0,1.005)

            ax.xaxis.offsetText.set_family(font["family"])
            ax.yaxis.offsetText.set_family(font["family"])

            # Set font of ticks
            for tick in ax.get_xticklabels():
                tick.set_fontname(font["family"])
                
            for tick in ax.get_yticklabels():
                tick.set_fontname(font["family"])

            # Save plot as pdf
            fig.savefig("../3d/benchmarks/benchmark-" + str(b+1) + "-" + str(s) + "/roofline_plot-" + str(b+1) + "-" + str(s) + ".pdf")
            plt.close(fig)
