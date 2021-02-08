#mpirun -np 2 ./RGFSolver 42 3 data/A_126_126_ns42_nt3.dat


ns=$1
nt=$2
no=$3

#folder_path=/home/gaedkem/spat_temp_m_pardiso/matrices_theta/theta_noise_5/theta_5_-10_2.5_1_ns${ns}_nt${nt}
#folder_path_data=/home/gaedkem/spat_temp_m_pardiso/data/temperature_data_2019/ns${ns}_nt${nt}

folder_path=/home/x_gaedkelb/spat_temp_m_pardiso/matrices_theta/theta_noise_5/theta_5_-10_2.5_1_ns${ns}_nt${nt}
folder_path_data=/home/x_gaedkelb/spat_temp_m_pardiso/data/temperature_data_2019/ns${ns}_nt${nt}

nb=2

echo "./RGFSolver ${folder_path} ${ns} ${nt} ${nb} ${folder_path_data} ${no} >${folder_path}/RGF_output_sel_inv.txt"
CUDA_VISIBLE_DEVICES="0" ./RGFSolver ${folder_path} ${ns} ${nt} ${nb} ${folder_path_data} ${no} >${folder_path}/RGF_output_sel_inv.txt
#mpirun ./RGFSolver ${folder_path} ${ns} ${nt} ${nb} ${folder_path_data} ${no} >${folder_path}/RGF_output_sel_inv.txt
