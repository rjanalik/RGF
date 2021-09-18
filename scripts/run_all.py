import glob
import os
import re
from tqdm import tqdm
import subprocess
OVERWRITE_RESULTS = True

def get_list_of_experiments(base_path):
    """
    returns list of (no, ns, nt) experiments
    """
    list_of_y_files = list()
    experiments = list()
    for (dirpath, dirnames, filenames) in os.walk(base_path):
        y_files = [file for file in filenames if file.startswith('y_')]
        # path_of_y_files += [os.path.join(dirpath, file) for file in filenames if file[0]=="y"]
        if y_files:
            y_file = y_files[0]
            pieces =  y_file.split("_")
            no = pieces[1]
            ns = int(re.sub('\D', '', pieces[3]))
            nt = int(re.sub('\D', '', pieces[4]))
            experiments.append((no, ns, nt))
    return experiments

if __name__ == '__main__':
    base_path = "../data/input/ghcn/2019/spatio_temporal/"
    experiments = get_list_of_experiments(base_path)
    print("Running %d experiments" % len(experiments))

    if(OVERWRITE_RESULTS):
        with open('../results/ghcn/results.csv', "w") as results_file:
            results_file.write('Date/Time,no,ns,nt,t_factorise,flops_factorize,t_solve,flops_solve,t_inv,flops_inv')
    for no, ns, nt in tqdm(experiments):
        print("Running no=%s, ns=%s nt=%s  experiments" % (no, ns, nt))
        error = subprocess.run(['./run_single.sh', 'no='+str(no), 'ns='+str(ns), 'nt='+str(nt)], shell=True, capture_output=True)
