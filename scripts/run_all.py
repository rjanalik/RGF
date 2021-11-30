import glob
import os
import re
from tqdm import tqdm
import subprocess
import argparse

parser = argparse.ArgumentParser(description='Run all available experiments in the data folder',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-o', '--overwrite', help='overwrite file', action='store_true')
args = parser.parse_args()
#OVERWRITE_RESULTS = args.overwrite
OVERWRITE_RESULTS = True

def get_list_of_experiments(base_path):
    """
    returns list of (no, ns, nt) experiments
    """
    list_of_y_files = list()
    experiments = list()
    for (dirpath, dirnames, filenames) in os.walk(base_path):
        # y_files = [file for file in filenames if file.startswith('rhs')]
        # # path_of_y_files += [os.path.join(dirpath, file) for file in filenames if file[0]=="y"]
        for folder in dirnames:
        # if y_files:
        #     y_file = y_files[0]
            pieces =  folder.split("_")
        #     no = pieces[1]
            ns = int(re.sub('\D', '', pieces[0]))
            nt = int(re.sub('\D', '', pieces[1]))
            nb = int(re.sub('\D', '', pieces[2]))
            print(nt)
            experiments.append((ns, nt, nb))
        return experiments

if __name__ == '__main__':
    base_path = "../data/input/tests/spatio_temporal/"
    experiments = get_list_of_experiments(base_path)
    print("Running %d experiments" % len(experiments))
    if(OVERWRITE_RESULTS):
        with open('../results/tests.csv', "w") as results_file:
            results_file.write('Date/Time\t\tns\tnt\tnb\tt_factorize\tflops_factorize\tt_solve\tflops_solve\tt_inv\tflops_inv\tRGF_Version\n')
    for ns, nt, nb in tqdm(experiments):
        print("Running ns=%s, nt=%s nb=%s  experiments" % (ns, nt, nb))
        result = subprocess.run(['bash ./run_single.sh '+str(ns)+' '+str(nt)+' '+str(nb)], shell=True, capture_output=True)
        print(result.stdout)
        print(result.stderr)
