ns=${1:-3}
nt=${2:-3}
nb=${3:-2}
set -x
bash ./C/generate_tests.sh ${ns} ${nt} ${nb}
bash ./run_single.sh ${ns} ${nt} ${nb}
