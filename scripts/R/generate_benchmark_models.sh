#!/usr/bin/env bash
set -euo pipefail

# for ns in 4002
# for ns in {1..10..1}
ns=20
# do
    for nt in {0..730..10}
    # for nt in {3}
    do
	if [[ "$nt" == 0 ]]
       	then
		# Sequential
		# Rscript generate_spatial_temporal_model.R -s $ns -t $nt 
		# Parallel
		Rscript generate_spatial_temporal_model.R -s $ns -t $(($nt+1)) &
	else
		# Sequential
		# Rscript generate_spatial_temporal_model.R -s $ns -t $nt 
		# Parallel
		Rscript generate_spatial_temporal_model.R -s $ns -t $nt &
	fi
    done
# done
