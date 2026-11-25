import os
import sys
from pathlib import Path
# Add parent directory to path so imports work regardless of where script is run from
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plotting.loading_functions import *
from plotting.plotting_functions import *
from models_main.helpers import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mc
from datetime import datetime, timedelta
import firedrake as df
from firedrake.pyplot import tripcolor
import cmcrameri.cm as cmc
import matplotlib.colors as colors
import fiona

fig_mode = ["main", "turb_vs_transition", "Re_omega", "three_catchments", "more_smoothing"][0]

# define data paths relative to this script's location
script_dir     = Path(__file__).resolve().parent
project_dir    = script_dir.parent
domain_dir     = project_dir / "domain"
outline_file   = domain_dir / "domain_outline.gpkg"
outline_small  = domain_dir / "domain_smaller.gpkg"
outline_Mank   = domain_dir / "domain_Mankoff_k09_5glaciers.gpkg"
data_dir       = project_dir /  "Greenland_data"
bedmachine_dir = data_dir / "bedmachine"
bed_file       = str(bedmachine_dir / "BedMachineGreenland-v5_bed_smooth_sig5.nc")
thick_file     = str(bedmachine_dir / "BedMachineGreenland-v5_thickness_smooth_sig5.nc")
fig_dir        = project_dir / "plotting" / "main_figures"

draw_Q = True
if fig_mode == "main":
    run_indices = [312,706]
    model_labels  = ["Baseline", "Reduced\nsheet flow"]
    fig_height = 13
    outline_paths = [outline_file, outline_file]
    figname = "f03"
elif fig_mode == "turb_vs_transition":
    run_indices = [602,706]
    model_labels  = ["Turbulent", "Reduced\nsheet flow"]
    fig_height = 13
    outline_paths = [outline_file, outline_file]
    figname = "fS12"
elif fig_mode == "Re_omega":
    run_indices = [706]
    model_labels  = ["Reduced\nsheet flow"]
    fig_height = 8
    outline_paths = [outline_file]
    figname = "fS13"
    draw_Q = False
elif fig_mode == "more_smoothing":
    run_indices = [718,706]
    model_labels  = ["More\nsmoothing", "Reduced\nsheet flow"]
    fig_height = 13
    outline_paths = [outline_file, outline_file]
    figname = "fS24"
elif fig_mode == "three_catchments":
    run_indices = [706,711,712]
    model_labels  = ["Red. sheet\nflow", "Smaller\ndomain", "Larger\ndomain"]
    fig_height = 19
    outline_paths = [outline_file, outline_small, outline_Mank]
    figname = "fS20"

# minimum discharge at which to draw channels
Q_min     = 15

# plotting parameters
plt.rcParams['font.size'] = 32
annotate_fs = 32

base_cmap   = cmc.lajolla_r
cmap_colors = base_cmap(np.linspace(0.1, 1.0, 256))
capped_map_Q = colors.ListedColormap(cmap_colors)

base_cmap = cmc.lapaz
cmap_colors = base_cmap(np.linspace(0.2, 0.85, 256))
capped_map_q = colors.ListedColormap(cmap_colors)
# "BrBG" for laminar vs turbulent, Re*omega
# capped_map_q = "BrBG"



def get_Q_lines(Q_vals, ax, cax, coords):
    # Get cell-vertex connectivity using CG function space
    V_CG = df.FunctionSpace(smesh_, "CG", 1)
    cells = V_CG.cell_node_list
    # Create line segments for mesh edges, each colored by Q value
    lines = []
    colors = []
    # Determine if cells are edges (2 vertices) or triangles (3 vertices)
    num_vertices_per_cell = cells.shape[1]
    for cell_idx, cell in enumerate(cells):
        q_val = Q_vals[cell_idx]
        vertices = cell
        if num_vertices_per_cell == 2:
            # Edges: just draw the single line segment
            p1 = coords[vertices[0]]
            p2 = coords[vertices[1]]
            lines.append([p1[:2], p2[:2]])
            colors.append(q_val)
        elif num_vertices_per_cell == 3:
            # Triangles: draw all 3 edges
            edges = [(vertices[0], vertices[1]), (vertices[1], vertices[2]), (vertices[2], vertices[0])]
            for edge in edges:
                p1 = coords[edge[0]]
                p2 = coords[edge[1]]
                lines.append([p1[:2], p2[:2]])
                colors.append(q_val)
    # Create LineCollection and add to plot
    # lc = mc.LineCollection(lines, cmap='Greys',
                        #    norm=plt.Normalize(vmin=np.nanmin(Q_vals)*(1), vmax=Q_min*1000))
    lc = mc.LineCollection(lines, cmap=capped_map_Q, norm=plt.Normalize(vmin=Q_min, vmax=140), joinstyle='round', capstyle='round')
                        #    norm=plt.Normalize(vmin=np.nanmin(Q_vals)*(1), vmax=Q_min))
    lc.set_array(np.array(colors))
    lc.set_linewidth(3)
    cbar = fig.colorbar(lc, cax=cax, label=r"$Q\,(\mathrm{m^3 s^{-1}})$")
    # manually define position and size of colorbar
    pos = cax.get_position()
    cax.set_position([
        pos.x0 + 0.2 * pos.width,
        pos.y0 + 0.25 * pos.height,
        0.2*pos.width,
        0.5*pos.height
    ])

    return lc

# Create figure and plot
Ly = np.zeros(len(run_indices))
for (i,r) in enumerate(run_indices):
    tp = str(project_dir / "parameter_runs" / f"run_{r}" / "time_series.h5")
    mesh_, smesh_ = get_meshes(tp)
    meshx, meshy = zip(*get_coordinates(mesh_, "CG", 1))
    Ly[i] = np.max(meshy) - np.min(meshy)

fig = plt.figure(figsize=(20, fig_height))
# height_ratios = Ly/Ly[0]
gs = fig.add_gridspec(len(run_indices), 3, left=0.15, right=0.91, top=0.9, bottom=0, hspace=-0.1, wspace=0, width_ratios=[1,1,0.12], height_ratios=Ly)

for (i_r,(run_index,outline_path)) in enumerate(zip(run_indices,outline_paths)):
    # load model output
    timeseries_path = str(project_dir / "parameter_runs" / f"run_{run_index}" / "time_series.h5")
    mesh_, smesh_ = get_meshes(timeseries_path)
    coords = smesh_.coordinates.dat.data  # mesh coordinates
    us_raw, m_raw, phi_raw, h_raw, q_raw, Q_raw, _, n_idx, _ = load_model_output(timeseries_path)
    for (i_s,season) in enumerate(["July 2020", "January 2021"]):
        idx = {'January 2021': int(7.041*365/2), 'July 2020': int(6.54*365/2)}[season]
        ax = fig.add_subplot(gs[i_r,i_s])

        # extract q for the specified timestep
        V   = df.FunctionSpace(mesh_, "CG", 1)
        q_vec_x = q_raw[::2, idx]
        q_vec_y = q_raw[1::2, idx]
        q = df.Function(V)
        if fig_mode == "Re_omega":
            q.dat.data[:] = np.sqrt(q_vec_x**2 + q_vec_y**2) / (3600*24*365) /1.787e-6 / 2000  # Re*omega
            cl = tripcolor(q, axes=ax, cmap=cmc.vik, norm=colors.LogNorm(vmin=1/50, vmax=50))
            clabel = r"Re $\omega$"
        else:
            q.dat.data[:] = np.sqrt(q_vec_x**2 + q_vec_y**2) / (3600*24*365)
            cl = tripcolor(q, axes=ax, cmap=capped_map_q, norm=colors.LogNorm(vmin=1e-5, vmax=0.1))
            clabel = r"$q\,(\mathrm{m^2 s^{-1}})$"

        # phi = df.Function(V)
        # phi.dat.data[:] = phi_raw[:,idx]
        # gradphi = df.max_value(df.dot(df.grad(phi),df.grad(phi)),1e-15)**(1/2)
        # keff = df.Function(V).interpolate(q / gradphi)


        # h = df.Function(V)
        # h.dat.data[:] = h_raw[:,idx]
        # h.dat.data[:] = 910*9.81*H.dat.data[:] - (phi_raw[:,idx] - 1000*9.81*B.dat.data[:])
        # h.interpolate(q / )

        # cl = tripcolor(h, axes=ax, cmap=cmc.lapaz, vmin=0, vmax=3e6) #, norm=colors.LogNorm(vmin=1e-3, vmax=3))
        # cl = tripcolor(keff, axes=ax, cmap=cmc.lapaz, norm=colors.LogNorm(vmin=1e-5, vmax=0.1))
        if (i_s == 0) & (i_r ==0):
            cax1 = fig.add_subplot(gs[0, 2])
            cbar = fig.colorbar(cl, cax=cax1, label=clabel) # cax= ... not ax=... !

            # manually define position and size of colorbar
            pos = cax1.get_position()
            cax1.set_position([
                pos.x0 + 0.2 * pos.width,
                pos.y0 + 0.25 * pos.height,
                0.2*pos.width,
                0.5*pos.height
            ])
        if (i_s == 0) & (i_r == len(run_indices)-1):
            # scale bar
            scale_length = 10e3
            x0 = -2.23e5
            y0 = -2.568e6
            line_x = [x0,x0+scale_length]
            line_y = [y0,y0]
            ax.plot(line_x, line_y, color="black", lw=7)
            ax.annotate(f"{int(scale_length*1e-3)} km", [x0-6e3,line_y[0]+5e3], c="black", fontsize=annotate_fs)

        # domain outline
        gdf = gpd.read_file(outline_path)
        gdf.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=2, label='Domain')

        if draw_Q:
            # extract Q for the specified timestep
            Q_vals = Q_raw[:, idx] / (3600*24*365)
            Q_vals[np.where(Q_vals<Q_min)] = np.nan
            cax2 = fig.add_subplot(gs[1, 2])
            lc = get_Q_lines(Q_vals, ax, cax2, coords)
            ax.add_collection(lc)

            meshx, meshy = zip(*get_coordinates(mesh_, "CG", 1))

        min_x, max_x = np.min(meshx), np.max(meshx)
        min_y, max_y = np.min(meshy), np.max(meshy)
        remove_axes(ax, ylims=(min_y,max_y), xlims=(min_x, min_x+0.5*(max_x-min_x)))

        r_bed = get_topo_raster(bed_file, outline_path)
        r_thick = get_topo_raster(thick_file, outline_path)
        r_surf  = r_bed+r_thick
        data    = r_surf.data.squeeze()
        x, y    = r_surf.coords(grid=False)
        cs = ax.contour(x, np.flip(y), data,
                        levels=np.arange(0, 2000, 100),
                        colors="black",
                        linewidths=0.7,
        )
        labels = ax.clabel(cs, fmt="%d m", fontsize=annotate_fs*0.8, levels=[1000,1400,], rightside_up=False)

        if i_s == 0:
            model_label = model_labels[i_r]
            if len(run_indices) > 1:
                x_offset = -0.52
            else:
                x_offset = -0.42
            ax.annotate(model_label, xy=(x_offset, 0.5), xycoords="axes fraction", annotation_clip=False, fontsize=annotate_fs, fontweight="bold", va="center")
        if i_r == 0:
            if len(run_indices) > 1:
                y_offset = 1.15
            else:
                y_offset = 1.07
            ax.annotate(season, xy=[0.6,y_offset], xycoords="axes fraction", annotation_clip=False, fontsize=annotate_fs, fontweight="bold", ha="center", va="center")
        ax.annotate(f"{idx_to_letter(i_r*2+i_s)}", xy=(0.03,1.05), xycoords="axes fraction", fontsize=annotate_fs*0.9, fontweight="bold")

plt.savefig(fig_dir / f"{figname}.jpg")
