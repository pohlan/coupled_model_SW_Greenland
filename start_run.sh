#!/bin/bash
# chmod u+x start_run.sh

target_directory='parameter_runs/'
run_index=706

# make a new directory and copy this file so that the input parameters are recorded at the start of the simulation
new_dir=$target_directory'run_'$run_index/
mkdir -p $new_dir
cp start_run.sh $new_dir'input_run'$run_index'.sh'

# parameters
e_v=0.0001
l_r=5
h_r=1.0
k_s=0.0006
k_c=0.1
l_c=10
beta2=3.7e5
p=1.2
q=1.0
transition=true
alpha_s=1.5
beta_s=1.5
omega=0.0005
As_factor=6
sig_topo=5
melt_input=MAR_monthly
m_basal=0.02
moulins=false
t_end=10
beta2_inversion=false
friction_heat=false

# record parameters in a txt file
cat > "$new_dir/parameters.txt" <<EOF
run_index=$run_index
e_v=$e_v
l_r=$l_r
h_r=$h_r
k_s=$k_s
k_c=$k_c
l_c=$l_c
beta2=$beta2
p=$p
q=$q
transition=$transition
alpha_s=$alpha_s
beta_s=$beta_s
omega=$omega
As_factor=$As_factor
sig_topo=$sig_topo
melt_input=$melt_input
m_basal=$m_basal
moulins=$moulins
t_end=$t_end
friction_heat=$friction_heat
EOF

# options that are true or false
options="--t_end $t_end --e_v $e_v --l_r $l_r --h_r $h_r --k_s $k_s --k_c $k_c --l_c $l_c --beta2 $beta2 --p $p --q $q --run_index $run_index --alpha_s $alpha_s --beta_s $beta_s --omega $omega --As_factor $As_factor --sig_topo $sig_topo --m_basal $m_basal --melt_input $melt_input"
if [ "$transition" = true ]; then
  options="$options --transition"
fi
if [ "$moulins" = true ]; then
  options="$options --moulins"
fi
if [ "$friction_heat" = true ]; then
  options="$options --friction_heat"
fi
if [ "$beta2_inversion" = true ]; then
  options="$options --beta2_inversion"
fi
# python -u inversion.py $options > parameter_runs/log_run_inversion_$run_index.txt
nohup python -u run_coupled_model.py $options > parameter_runs/log_run_$run_index.txt &
