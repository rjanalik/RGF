#mpirun -np 2 ./RGFSolver 42 3 data/A_126_126_ns42_nt3.dat


: ${ns:=1212}
: ${nt:=20}
: ${no:=236680}
# ns=$1
# nt=$3
# no=$3

#folder_path=/home/x_gaedkelb/RGF/data/ns${ns}
folder_path=/home/x_pollakgr/RGF/data/input/ghcn/2019/spatio_temporal/ns${ns}_nt${nt}

nb=2
year=2019

echo "GV100"
echo "GPU 0 ./main ${folder_path} ${ns} ${nt} ${nb} ${no} >${folder_path}/RGF_output_sel_inv.txt"

#CUDA_VISIBLE_DEVICES="0" ./mainConstInd ${folder_path} ${ns} ${nt} ${nb} ${folder_path_data} ${no} >${folder_path}/RGF_output.txt
CUDA_VISIBLE_DEVICES="0" /home/x_pollakgr/RGF/build/bin/main ${folder_path} ${ns} ${nt} ${nb} ${folder_path_data} ${no}
#mv ${folder_path_data}/log_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_year${year}.dat ${folder_path_data}/log_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_year${year}_gpu0.dat