#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH -c 4
#SBATCH --partition cpu
#SBATCH --cpus-per-node=2
#SBATCH --mem=64G
#SBATCH --job-name=dedup_th_tagged_1
#SBATCH -o log_1.out

NAME=$1
SM_FILE=$2
IDX=$3
now=$(date +"%T")
echo "Current time : $now"
py-spy record -o ${NAME}_${IDX}.svg --rate 1 -- python /hpc/compgen/projects/gw_cfdna/gw_cyclomics/analysis/lchen/FragmentOmics/cyclomicsConsensus/dedup_tidehunter_tagged.py --read_bam /hpc/compgen/projects/gw_cfdna/gw_cyclomics/analysis/lchen/FragmentOmics-output/dedup_by_coordinate/${SM_FILE} --out_path /hpc/compgen/projects/gw_cfdna/gw_cyclomics/analysis/lchen/FragmentOmics-output/dedup_by_coordinate/ --prefix ${NAME}_${IDX} --ref /hpc/cog_bioinf/GENOMES.old/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta
#python /hpc/compgen/projects/gw_cfdna/gw_cyclomics/analysis/lchen/FragmentOmics/dedup.py --read_bam /hpc/compgen/projects/gw_cfdna/gw_cyclomics/analysis/lchen/concall_all/concall_github/concall-working/concall/output/MAR6557_b21_99L/05_aggregated/MAR6557_b21_99L_tide_rotated.tagged.sorted.bam --SM MAR6557_b21_99L --out_path /hpc/compgen/projects/gw_cfdna/gw_cyclomics/analysis/lchen/FragmentOmics/data --prefix MAR6557_b21_99L --ref /hpc/cog_bioinf/GENOMES.old/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta
dmesg > /hpc/compgen/projects/gw_cfdna/gw_cyclomics/analysis/lchen/FragmentOmics-output/dedup_by_coordinate/dmesg_log/dmesg-${NAME}_${IDX}_$!.log
now=$(date +"%T")
echo "Current time : $now"
