#mpirun -np 2 ./RGFSolver 42 3 data/A_126_126_ns42_nt3.dat

ns=42
nt=3
no=128
# nu = ns * nt = 42*3 = 126

#folder_path=/home/x_gaedkelb/RGF/data/ns${ns}
folder_path=/home/x_grpollak/RGF/data/input/ghcn/2019/spatio_temporal/ns${ns}_nt${nt}

nb=2
year=2019

echo "GV100"
# echo "GPU 0 ./mainConstInd ${folder_path} ${ns} ${nt} ${nb} ${no} >${folder_path}/RGF_output_sel_inv.txt"

CUDA_VISIBLE_DEVICES="0" ~/RGF/build/bin/main --path ${folder_path} --ns ${ns} --nt ${nt} --nb ${nb} --no ${no}
#mv ${folder_path}/L_factor_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}.dat ${folder_path}/L_factor_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_gpu0.dat
#mv ${folder_path_data}/x_sol_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_year${year}.dat ${folder_path_data}/x_sol_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_year${year}_gpu0.dat
#mv ${folder_path_data}/log_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_year${year}.dat ${folder_path_data}/log_RGF_ns${ns}_nt${nt}_nb${nb}_no${no}_year${year}_gpu0.dat

