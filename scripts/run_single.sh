#!/bin/bash
#mpirun -np 2 ./RGFSolver 42 3 data/A_126_126_ns42_nt3.dat
source /opt/intel/oneapi/mkl/latest/env/vars.sh  intel64
export LD_LIBRARY_PATH=/home/x_pollakgr/RGF/applications/magma-2.5.4/lib:$LD_LIBRARY_PATH
export MAGMA_DIR=/home/x_pollakgr/RGF/applications

ns=${1:-3}
nt=${2:-4}
nb=${3:-2}
year=${5:-2019}

while [ $# -gt 0 ]; do
   if [[ $1 == *"b'--"* ]]; then
        param="${1/--/}"
        declare $param="$1"
        echo $1 $2 "success"
   fi
  shift
done

#folder_path=/home/x_gaedkelb/RGF/data/ns${ns}
if [[ "$nt" == 1 ]];
then
	folder_path=/home/x_pollakgr/RGF/data/input/tests/spatial/ns${ns}_nt${nt}_nb${nb}
else
	folder_path=/home/x_pollakgr/RGF/data/input/tests/spatio_temporal/ns${ns}_nt${nt}_nb${nb}
fi

echo "GV100"
#echo "GPU 1 ./main -path ${folder_path} -ns ${ns} -nt ${nt} -nb ${nb} -no ${no} >${folder_path}/RGF_output_sel_inv.txt"
#set -x
export CUDA_VISIBLE_DEVICES=0
echo "CUDA_VISIBLE_DEVICES = " $CUDA_VISIBLE_DEVICES
CUDA_VISIBLE_DEVICES="1" /home/x_pollakgr/RGF/build/bin/main --path ${folder_path} --ns ${ns} --nt ${nt} --nb ${nb}
# LD_LIBRARY_PATH=/home/x_pollakgr/RGF/external/magma-2.5.4/lib:$LD_LIBRARY_PATH
#/home/x_pollakgr/RGF/build/bin/main 3 3 0 "/home/x_pollakgr/RGF/data/A_9_9_ns3_nt3.dat" "/home/x_pollakgr/RGF/data/rhs9.txt"
#/home/x_pollakgr/RGF/build/bin/main 3 3 2 "/home/x_pollakgr/RGF/data/A_11_11_ns3_nt3_nd2.dat" "/home/x_pollakgr/RGF/data/rhs11.txt"
#/home/x_pollakgr/RGF/build/bin/main 4 4 2 "/home/x_pollakgr/RGF/data/input/tests/spatio_temporal/ns_4_nt4/ns_4_nt4_nb2.mat" "/home/x_pollakgr/RGF/data/rhs18.txt"
#gdb --args /home/x_pollakgr/RGF/build/bin/main 10 10 5 "/home/x_pollakgr/RGF/data/matrixRGF.mat" "/home/x_pollakgr/RGF/data/rhs105.txt"
#/home/x_pollakgr/RGF/build/bin/main 10 10 5 "/home/x_pollakgr/RGF/data/matrixRGF.mat" "/home/x_pollakgr/RGF/data/rhs105.txt"
#/home/x_pollakgr/RGF/build/bin/main 3 3 2 "/home/x_pollakgr/RGF/data/matrixRGF_ns3_nt3_nd2.mat" "/home/x_pollakgr/RGF/data/rhs11.txt"
#/home/x_pollakgr/RGF/build/bin/main 1002 16 2 "/home/x_pollakgr/RGF/data/input/ghcn/2019/spatio_temporal/ns1002_nt16/Qxy_R_16034.dat" "/home/x_pollakgr/RGF/data/input/ghcn/2019/spatio_temporal/ns1002_nt16/bxy_R_16034_1.dat"

#gdb --args /home/x_pollakgr/RGF/build/bin/main --path /home/x_pollakgr/RGF/data/input/ghcn/2019/spatio_temporal/ns4002_nt10 --ns 4002 --nt 10 --nb 2 --no 118459
#/home/x_pollakgr/RGF/build/bin/main --path /home/x_pollakgr/RGF/data/input/ghcn/2019/spatio_temporal/ns4002_nt100 --ns 4002 --nt 100 --nb 2 --no 1188566
#/home/x_pollakgr/RGF/build/bin/main --path /home/x_pollakgr/RGF/data/input/ghcn/2019/spatio_temporal/ns1002_nt16 --ns 1002 --nt 16 --nb 2 --no 188242
#nvprof -f -o result.nvvp /home/x_pollakgr/RGF/build/bin/main --path /home/x_pollakgr/RGF/data/input/ghcn/2019/spatio_temporal/ns4002_nt10 --ns 4002 --nt 10 --nb 2 --no 118459
#mv ${folder_path_data}/log_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_year${year}.dat ${folder_path_data}/log_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_year${year}_gpu0.dat
