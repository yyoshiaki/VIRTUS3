#!/bin/bash                                       
#SBATCH --ntasks=1
#SBATCH --partition=day
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20GB
#SBATCH --time=10:00:00
#SBATCH --mail-user=heng-le.chen@yale.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/gpfs/gibbs/pi/hafler/hc865/VIRTUS3/slogs/test_virtus_%j.out

set -x
module purge
module load miniconda
conda activate /vast/palmer/pi/hafler/hc865/conda_envs/virtus3/
module load CellRanger/9.0.1
module load Salmon/1.4.0-gompi-2020b
module load SAMtools/1.21-GCC-12.2.0

python /gpfs/gibbs/pi/hafler/hc865/VIRTUS3/src/virtus3.py \
    --fastqs /home/hc865/palmer_scratch/SRR16976513_data \
    --chemistry_cr ARC-v1 \
    --sample SRR16976513 \
    --lib_alevin="-l ISR --chromiumV3" \
    --output /gpfs/gibbs/pi/hafler/hc865/VIRTUS3/test_output \
    --index_human /gpfs/gibbs/pi/hafler/hc865/VIRTUS3/data/refdata-gex-GRCh38-2024-A \
    --index_virus /gpfs/gibbs/pi/hafler/hc865/VIRTUS3/data/NC_007605.1_CDS_EBER12_salmon_index \
    --tgMap /gpfs/gibbs/pi/hafler/hc865/VIRTUS3/data/NC_007605.1_CDS_EBER12.tgMap.tsv \
    --cellranger /vast/palmer/apps/avx2/software/CellRanger/9.0.1/cellranger \
    --salmon /vast/palmer/apps/avx2/software/Salmon/1.4.0-gompi-2020b/bin/salmon \
    --cores 40