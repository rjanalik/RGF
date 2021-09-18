#!/bin/bash
#mpirun -np 2 ./RGFSolver 42 3 data/A_126_126_ns42_nt3.dat


no=$1
ns=$2
nt=$3
# environment=${environment:-production}
# school=${school:-is out}

no=${1:-236680}
ns=${2:-1212}
nt=${3:-20}
nb=${4:-2}
year=${5:-2019}

while [ $# -gt 0 ]; do
   if [[ $1 == *"b'--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        echo $1 $2 "success"
   fi
  shift
done

#folder_path=/home/x_gaedkelb/RGF/data/ns${ns}
folder_path=/home/x_pollakgr/RGF/data/input/ghcn/${year}/spatio_temporal/ns${ns}_nt${nt}

echo "GV100"
echo "GPU 0 ./main ${folder_path} ${ns} ${nt} ${nb} ${no} >${folder_path}/RGF_output_sel_inv.txt"
CUDA_VISIBLE_DEVICES="0" /home/x_pollakgr/RGF/build/bin/main ${folder_path} ${ns} ${nt} ${nb} ${folder_path_data} ${no}
#mv ${folder_path_data}/log_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_year${year}.dat ${folder_path_data}/log_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_year${year}_gpu0.dat
