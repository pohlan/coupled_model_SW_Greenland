# coupled_model_SW_Greenland

This code was adapted from original code by Doug Brinkerhoff, written for the paper [Brinkerhoff et al., 2021](https://doi.org/10.1017/jog.2020.112).

## Software requirements

The major requirement of this repository is firedrake. Detailed instructions for installation are provided [on their website](https://www.firedrakeproject.org/install.html).

For other python libraries see `requirements.txt`.

## Data requirements

- Bed and surface topography as a netcdf file
- Melt input: hard-coded to use MAR monthly or daily runoff; daily data can be downloaded from [here](https://doi.org/10.5281/zenodo.19691263), monthly data used to be available [here](http://ftp.climato.be/fettweis/) (currently not?)
- For moulin runs: Surface catchments, available for download in the supporting information of [Yang and Smith (2016)](https://doi.org/10.1002/2016JF003927)


[Note: this is not fully streamlined yet]

## Running the model

All the 'meat' of the model is in `models_main/coupled_model.py`. It can be used independently and does not depend on the type of input data, mesh or other setups that are built around it in the other files. The actual time stepping is done in `run_coupled_model.py`.

To use the model with a different geometry and/or input data, it is easiest to take `run_coupled_model.py` as an example and adapt it accordingly.

To start a model run with this setup here, the easiest is to open `start_run.sh`, modify parameters as desired, choose an identifier (`run_index`) and then execute it with `./start_run.sh`. It will call `run_coupled_model.py` and save all results in `parameter_runs/run_{run_index}/`. Because the model runs take several hours, the python script is called via nohup so that the process is not killed as easily. Printing of timesteps is logged in `parameter_runs/log_run_{run_index}.txt`.

In order to start several runs in one go, e.g. for a sensitivity analysis, use any of the `start_sensitivity_runs_*.sh` scripts. Starting more than 6-8 runs at once is not recommended, it will likely max out the memory for a 64 GB RAM computer.
