#!/usr/bin/env bash
set -euo pipefail

# for ns in 4002
# for ns in {1..10..1}
ns=20
# do
    for nt in {1..20..1}
    # for nt in {3}
    do
        Rscript generate_spatial_temporal_model.R -s $ns -t $nt 
    done
# done
