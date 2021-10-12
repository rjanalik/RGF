#!/bin/bash
#mpirun -np 2 ./RGFSolver 42 3 data/A_126_126_ns42_nt3.dat


no=$1
ns=$2
nt=$3
# environment=${environment:-production}
# school=${school:-is out}

no=${1:-11790}
ns=${2:-1002}
nt=${3:-1}
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
if [[ "$nt" == 1 ]]; 
then
	folder_path=/home/x_pollakgr/RGF/data/input/ghcn/${year}/spatial/ns${ns}_nt${nt}
else
	folder_path=/home/x_pollakgr/RGF/data/input/ghcn/${year}/spatio_temporal/ns${ns}_nt${nt}
fi

echo "GV100"
#echo "GPU 1 ./main -path ${folder_path} -ns ${ns} -nt ${nt} -nb ${nb} -no ${no} >${folder_path}/RGF_output_sel_inv.txt"
set -x
# CUDA_VISIBLE_DEVICES="1" /home/x_pollakgr/RGF/build/bin/main --path ${folder_path} --ns ${ns} --nt ${nt} --nb ${nb} --no ${no}
export CUDA_VISIBLE_DEVICES="0"
# LD_LIBRARY_PATH=/home/x_pollakgr/RGF/external/magma-2.5.4/lib:$LD_LIBRARY_PATH
#/home/x_pollakgr/RGF/build/bin/main --path /home/x_pollakgr/RGF/data/input/ghcn/2019/spatio_temporal/ns4002_nt10 --ns 4002 --nt 10 --nb 2 --no 118459
nvprof -f -o result.nvvp /home/x_pollakgr/RGF/build/bin/main --path /home/x_pollakgr/RGF/data/input/ghcn/2019/spatio_temporal/ns4002_nt10 --ns 4002 --nt 10 --nb 2 --no 118459
#mv ${folder_path_data}/log_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_year${year}.dat ${folder_path_data}/log_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_year${year}_gpu0.dat
