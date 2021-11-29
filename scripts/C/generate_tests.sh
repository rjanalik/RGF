#!/usr/bin/env bash
set -euo pipefail

#!/usr/bin/env bash
set -euo pipefail

# for ns in 4002
# for ns in {1..10..1}
# ns = 500, nt = {50, 100, 200, 500}, nb = 5
# ns = 1000, nt = {50, 100, 200, 500}, nb = 5

HOME=/home/x_pollakgr/RGF
DATA_FOLDER=${HOME}/data/input/tests

ns=3
nb=2
# set -x
for nt in {0..4..1}
# for nt in {3}
do
let rows=$((${ns}*${nt}+${nb}))
    if [[ "$nt" == 0 ]]
        then
        # Sequential
        #/home/x_pollakgr/RGF/scripts/C/.inverse_laplace_dense -nx $ns -ny $nt -nz 1 -fname ns_${ns}_nt${nt}_nb${nb}.mat
        # Parallel
        OUTPUT_FOLDER=${DATA_FOLDER}/spatial/ns${ns}_nt${nt}_nb${nb}
        mkdir -p ${OUTPUT_FOLDER};
        ${HOME}/scripts/C/inverse_laplace_dense -nx=${ns} -ny=$(($nt+1)) -nz=1 -nb=${nb} -fname=${OUTPUT_FOLDER}/ns${ns}_nt${nt}_nb${nb}.mat &
    else
        # Sequential
        #/home/x_pollakgr/RGF/scripts/C/.inverse_laplace_dense -nx $ns -ny $nt -nz 1 -fname ns_${ns}_nt${nt}_nb${nb}.mat
        # Parallel
        OUTPUT_FOLDER=${DATA_FOLDER}/spatio_temporal/ns${ns}_nt${nt}_nb${nb}
        mkdir -p ${OUTPUT_FOLDER};
        ${HOME}/scripts/C/inverse_laplace_dense -nx=${ns} -ny=${nt} -nz=1 -nb=${nb} -fname=${OUTPUT_FOLDER}/ns${ns}_nt${nt}_nb${nb}.mat &
    fi
    # rhs files
    for i in {1..${rows}}
    do
        echo "${i}\n" > "${OUTPUT_FOLDER}/rhs${rows}.txt"
    done
done
