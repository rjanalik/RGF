#!/bin/bash
#mpirun -np 2 ./RGFSolver 42 3 data/A_126_126_ns42_nt3.dat
source /opt/intel/oneapi/mkl/latest/env/vars.sh  intel64
export LD_LIBRARY_PATH=/home/x_pollakgr/RGF/applications/magma-2.5.4/lib:$LD_LIBRARY_PATH
export MAGMA_DIR=/home/x_pollakgr/RGF/applications
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# ns = 500, nt = {50, 100, 200, 500}, nb = 5
# ns = 1000, nt = {50, 100, 200, 500}, nb = 5

ns=${1:-4}
nt=${2:-8}
nb=${3:-2}
nvvp_file=${4}
#year=${5:-2019}

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
	folder_path=/home/x_pollakgr/RGF/data/input/tests/spatial/ns${ns}_nt${nt}_nb${nb}
 else
	folder_path=/home/x_pollakgr/RGF/data/input/tests/spatio_temporal/ns${ns}_nt${nt}_nb${nb}
 fi

#echo "GPU 1 ./main -path ${folder_path} -ns ${ns} -nt ${nt} -nb ${nb} -no ${no} >${folder_path}/RGF_output_sel_inv.txt"
# CUDA_VISIBLE_DEVICES="1" /home/x_pollakgr/RGF/build/bin/main --path ${folder_path} --ns ${ns} --nt ${nt} --nb ${nb} --no ${no}
export CUDA_VISIBLE_DEVICES=0
echo "CUDA_VISIBLE_DEVICES = " $CUDA_VISIBLE_DEVICES
set -x
#CUDA_VISIBLE_DEVICES="0" /home/x_pollakgr/RGF/build/bin/main --path ${folder_path} --ns ${ns} --nt ${nt} --nb ${nb}
if [ -z "$nvvp_file" ]
then
    CUDA_VISIBLE_DEVICES="0" /home/x_pollakgr/RGF/build/bin/main --path ${folder_path} --ns ${ns} --nt ${nt} --nb ${nb} --results ${SCRIPT_DIR}/../results/tests.csv
else
    CUDA_VISIBLE_DEVICES="0" nvprof -f -o ${file}.nvvp /home/x_pollakgr/RGF/build/bin/main --path ${folder_path} --ns ${ns} --nt ${nt} --nb ${nb} --results ${SCRIPT_DIR}/../results/tests.csv
fi
# CUDA_VISIBLE_DEVICES="1" gdb --args /home/x_pollakgr/RGF/build/bin/main --path ${folder_path} --ns ${ns} --nt ${nt} --nb ${nb}
 # LD_LIBRARY_PATH=/home/x_pollakgr/RGF/external/magma-2.5.4/lib:$LD_LIBRARY_PATH
#/home/x_pollakgr/RGF/build/bin/main 3 3 0 "/home/x_pollakgr/RGF/data/input/Radim/A_9_9_ns3_nt3.dat" "/home/x_pollakgr/RGF/data/input/Radim/rhs9.txt"
#/home/x_pollakgr/RGF/build/bin/main 3 3 2 "/home/x_pollakgr/RGF/data/input/Radim/A_11_11_ns3_nt3_nd2.dat" "/home/x_pollakgr/RGF/data/input/Radim/rhs11.txt"
#/home/x_pollakgr/RGF/build/bin/main 4 2 2 "/home/x_pollakgr/RGF/data/input/tests/spatio_temporal/ns4_nt2_nb2/ns4_nt2_nb2.mat" "/home/x_pollakgr/RGF/data/input/tests/spatio_temporal/ns4_nt2_nb2/rhs10.txt"
#/home/x_pollakgr/RGF/data/input/tests/spatio_temporal/ns3_nt4_nb2/ns3_nt4_nb2.mat /home/x_pollakgr/RGF/data/input/tests/spatio_temporal/ns3_nt4_nb2/rhs14.txt
 #gdb --args /home/x_pollakgr/RGF/build/bin/main 10 10 5 "/home/x_pollakgr/RGF/data/matrixRGF.mat" "/home/x_pollakgr/RGF/data/rhs105.txt"
 #/home/x_pollakgr/RGF/build/bin/main 10 10 5 "/home/x_pollakgr/RGF/data/matrixRGF.mat" "/home/x_pollakgr/RGF/data/rhs105.txt"
 #/home/x_pollakgr/RGF/build/bin/main 3 3 2 "/home/x_pollakgr/RGF/data/matrixRGF_ns3_nt3_nd2.mat" "/home/x_pollakgr/RGF/data/rhs11.txt"4
#/home/x_pollakgr/RGF/build/bin/main --path /home/x_pollakgr/RGF/data/input/ghcn/2019/spatio_temporal/ns4002_nt100 --ns 4002 --nt 100 --nb 2 --no 1188566
#/home/x_pollakgr/RGF/build/bin/main --path /home/x_pollakgr/RGF/data/input/ghcn/2019/spatio_temporal/ns1002_nt16 --ns 1002 --nt 16 --nb 2 --no 188242
#nvprof -f -o result.nvvp /home/x_pollakgr/RGF/build/bin/main --path /home/x_pollakgr/RGF/data/input/ghcn/2019/spatio_temporal/ns4002_nt10 --ns 4002 --nt 10 --nb 2 --no 118459
#mv ${folder_path_data}/log_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_year${year}.dat ${folder_path_data}/log_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_year${year}_gpu0.dat

