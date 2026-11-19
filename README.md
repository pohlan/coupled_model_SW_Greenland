# coupled_model_SW_Greenland


## Software requirements

The major requirement of this repository is firedrake. Detailed instructions for installation are provided [on their website](https://www.firedrakeproject.org/install.html).

For other python libraries see `requirements.txt`.

## Data requirements

- Bed and surface topography as a netcdf file
- Melt input: hard-coded to use MAR monthly or daily runoff; daily data can be downloaded from [here](https://doi.org/10.5281/zenodo.19691263), monthly data used to be available [here](http://ftp.climato.be/fettweis/) (currently not?)
- For moulin runs: Surface catchments, available for download in the supporting information of [Yang and Smith (2016)](https://doi.org/10.1002/2016JF003927)


[Note: this is not fully streamlined yet]

## Running the model

To initiate a model run, the easiest is to open `start_run.sh`, modify parameters as desired, choose an identifier (`run_index`) and then execute it with `./start_run.sh`. It will call `run_coupled_model.py` and save all results in `parameter_runs/run_{run_index}/`. Because the model runs take several hours, the python script is called via nohup so that the process is not killed as easily. Printing of timesteps is logged in `parameter_runs/log_run_{run_index}.txt`.

In order to start several runs as for the sensitivity analysis, use any of the `start_sensitivity_runs_*.sh` scripts.

