#!/usr/bin/env bash
set -euo pipefail

for ns in {1..10..1}
do
    for nt in 20
    do
        Rscript generate_spatial_temporal_model.R -s $ns -t $nt &
    done
done
