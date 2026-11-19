import os
import sys
from pathlib import Path
# Add parent directory to path so imports work regardless of where script is run from
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plotting.loading_functions import *
from plotting.plotting_functions import *
from models_main.catchments import get_catchments
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import cmcrameri.cm as cmc
from firedrake.pyplot import tripcolor
from matplotlib.colors import LogNorm


# plotting parameters
color_map = cmc.managua
plt.rcParams['font.size'] = 29
annotate_fs = 29
lw = 2.5
clims = (20,300)

# define data paths relative to this script's location
script_dir      = Path(__file__).resolve().parent
project_path    = script_dir.parent
data_dir        = project_path / "Greenland_data"
bedmachine_dir  = data_dir
outline_path    = data_dir / "russel" / "russel_domain.gpkg"
catchment_dir   = data_dir / "Yang_Smith_2016" / "ds01"
catchments_file = catchment_dir / "catchments_LC80070132013231LGN00_B8.shp"
fig_path        = project_path / "plotting" / "main_figures"

year = 2020
xstart, xend = datetime(year,4,2), datetime(year,10,30)

# prepare figure, gridspecs
fig, ax = plt.subplots(figsize=(15,6))
fig.subplots_adjust(left=0.0, right=1.0, bottom=0.01, top=0.99)

# load model stuff not run-specific
run_index = 713
timeseries_path = str(project_path / "parameter_runs" / f"run_{run_index}" / "time_series.h5")
mesh_, smesh_ = get_meshes(timeseries_path)
df_moul, facet_functions, distributed_melt_mask = get_catchments(mesh_, bedmachine_dir, catchment_dir, data_dir, outline_path)

_, m_raw, _, _, _, _, _, n_idx, _ = load_model_output(timeseries_path)
with h5py.File(timeseries_path, 'r') as h5file:
    Qm_raw = h5file['topologies/firedrake_default_topology/dms/firedrake_dm_1_0_0_False_1/vecs/Qm/Qm'][()].T
dates_model, i_model = get_model_dates(n_idx, xstart, xend)
V     = df.FunctionSpace(mesh_, "CG", 1)
m     = df.Function(V)
m_distributed   = df.Function(V)
for (i_time, i_m) in enumerate(i_model):
    m_distributed.dat.data[:] = (m_raw[:,i_m]-0.02)/365*2 # remove basal melt 0.02 m/yr
    m_distributed.interpolate(m_distributed*distributed_melt_mask)
    m.dat.data[:] += Qm_raw[:, i_m]/365*2 + m_distributed.dat.data[:]

m.dat.data[m.dat.data_ro<1e-10] = np.nan
im = tripcolor(m, axes=ax, norm=LogNorm(vmin=1e-3, vmax=30), cmap=cmc.batlowK) # cmap=inferno
cbar = fig.colorbar(im, ax=ax, shrink=0.7)
cbar.set_label(r"$\dot{m}_\mathrm{surface}$" +r" ($\mathrm{m}$)" + "\nintegrated over 2020")

# format
ax.set_aspect('equal')
ax.set_xticks([])
ax.set_yticks([])
ax.spines[['right', 'left', 'top', 'bottom']].set_visible(False)
ax.set_ylim(-2.585e6,-2.47e6)
ax.set_xlim(-2.4e5,-2e4)

# domain outline
gdf = gpd.read_file(outline_path)
gdf.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=1, label='Domain')

# catchments
catchments = gpd.read_file(catchments_file).geometry
catchments = catchments.to_crs(gdf.crs)
catchments.plot(ax=ax, facecolor='none', edgecolor='grey', linewidth=0.3) #, alpha=0.5)
# get moulin coordinates that are exactly on the nodes
coords = hlp.get_coordinates(mesh_, "CG", 1)
meshx = coords[:,0]
meshy = coords[:,1]
x_moulin = np.zeros(len(df_moul.x))
y_moulin = np.zeros(len(df_moul.y))
id_mesh = np.zeros(len(df_moul.y), dtype=int)
i  = 0
for (mx,my) in zip(df_moul.x,df_moul.y):
    ii = np.argmin(np.sqrt((meshx-mx)**2+(meshy-my)**2))
    x_moulin[i] = meshx[ii]
    y_moulin[i] = meshy[ii]
    id_mesh[i] = ii
    i += 1
ax.scatter(x_moulin,y_moulin,2)

plt.savefig(fig_path / "fS15.jpg")
