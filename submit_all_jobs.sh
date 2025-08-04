#!/bin/bash

echo "=== Submitting Set 1: weighted_by_area with different nbins ==="
for sample in P2_CRC P5_CRC; do
  for nbins in 2 4 6 8 10; do
    echo "Submitting: $sample with weighted_by_area and nbins $nbins"
    sbatch /gpfs/commons/groups/innovation/jjoseph/enact-pipeline/run_individual_sample_batch.sh $sample weighted_by_area nbins$nbins
  done
done

echo ""
echo "=== Submitting Set 2: All methods with default nbins ==="
for sample in P2_CRC P5_CRC; do
  for method in naive weighted_by_area weighted_by_gene weighted_by_cluster; do
    echo "Submitting: $sample with $method and default nbins"
    sbatch /gpfs/commons/groups/innovation/jjoseph/enact-pipeline/run_individual_sample_batch.sh $sample $method default
  done
done

echo ""
echo "Total jobs submitted: 18 (10 Set 1 + 8 Set 2)"
