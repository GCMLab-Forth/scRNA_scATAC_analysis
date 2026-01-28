#!/bin/bash

# This is a shared resource so please be cautious and responsible!

# 230208 edialynas v0.1 sbatch skeleton script

#SBATCH --job-name=50x_scenic             # name of script, this is just for reporting/accounting purposes
#SBATCH --output=./50x_scenic.out          # standard output file
#SBATCH --error=./50x_scenic.err           # standard error file
#SBATCH --nodes=1                       # number of nodes to allocate, if your application does not run in parallel (MPI) mode set this to 1
#SBATCH --ntasks=40                    # number of cores to allocate, our nodes have 48 cores
#SBATCH --time=200:00:00                 # set a lim on the total run time, hrs:min:sec
#SBATCH --mem=70G                       # memory to allocate, our nodes have 256G, so set this to no more than 250G

# Remove the following two lines if no email notification is required
#SBATCH --mail-type=end,fail            # events to send mail on
#SBATCH --mail-user=marianna_stagaki@imbb.forth.gr    # mail recipient


WORKDIR=./                              # set this to your working directory
echo $WORKDIR
cd $WORKDIR
pwd;hostname;date                       # print-and log-working directory, hostname, start time & date, useful for job tracking and debugging


#!/bin/bash

# Directories and input files
TF_LIST="mm_mgi_tfs.txt"
MOTIF_DB="mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
ANNOTATIONS="motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
N_WORKERS=40

# Loop over all cluster matrices
for matrix in $(ls expr_matrices_per_condition/expression_matrix_cluster_Th1_Th17WT.csv); do   
    base=$(basename "$matrix" .csv)
    echo "🧬 Running 10 SCENIC iterations for: $base"
    
    # Run SCENIC 50 times with different seeds
    for i in {1..50}; do
        seed=$((42 + i))  # Change seed each iteration

        echo "  🔁 Iteration $i (seed $seed)..."

        # Step 1: GRN inference
        pyscenic grn --num_workers $N_WORKERS -o grn_${base}_rep${i}.tsv "$matrix" "$TF_LIST" --seed $seed

        # Step 2: Motif enrichment
        pyscenic ctx grn_${base}_rep${i}.tsv "$MOTIF_DB" \
            --annotations_fname "$ANNOTATIONS" \
            --expression_mtx_fname "$matrix" \
            --auc_threshold 0.15 \
            --mode "custom_multiprocessing" \
            --output regulons_${base}_rep${i}.csv \
            --num_workers $N_WORKERS

        # Step 3: AUCell scoring
        pyscenic aucell "$matrix" regulons_${base}_rep${i}.csv \
            -o auc_${base}_rep${i}.csv --num_workers $N_WORKERS

        echo "    ✅ Done rep $i for $base"
    done
done

echo "🎉 All SCENIC runs complete!"


