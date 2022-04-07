# run_script for main


ns=$1
#nt=3
nb=$2
nrhs=$3

#matrix_file="/home/x_gaedkelb/RGF_GPU_solve/RGF/data/A_9_9_ns3_nt3.dat"
matrix_file="/home/x_gaedkelb/RGF_GPU_solve/RGF/data/A_126_126_ns42_nt3.dat"
rhs_file="/home/x_gaedkelb/RGF_GPU_solve/RGF/data/b_9.dat"

#./main ${ns} ${nt} ${nd} ${matrix_file} 100
CUDA_VISIBLE_DEVICES="0" ./main ${ns} ${nb} ${nrhs}
