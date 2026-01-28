
from collections import defaultdict
import matplotlib.pyplot as plt 
import numpy as np
import csv
import ast
from ast import literal_eval
from collections import Counter
import pandas as pd 
from itertools import chain
# Initialize dictionaries to store AUCs, NES, and appearances

def process_regulon_file(filename, number_of_repetitions, condition, cluster, percentage = 0.3):
    # Process each run
    print(filename)
    tfs_aucs = defaultdict(list)
    tfs_nes = defaultdict(list)
    tfs_target_genes = defaultdict(list)
    for i in range(1, number_of_repetitions+1):
        with open(filename+str(i)+".csv", 'r') as f:
            lines = f.readlines()
            tf_max_nes = {}  # Temporary dictionary to store max NES per TF for this run
            for line in lines[3:]:
                tf_name = line.split(",")[0]
                auc = float(line.split(",")[2])
                nes = float(line.split(",")[3])
                target_genes_list = ast.literal_eval(" ".join(line.split("\"")[-2:-1]))
                target_genes = [gene[0] for gene in target_genes_list]
                # Update max NES for the TF if it's higher than the current max
                if tf_name not in tf_max_nes or nes > tf_max_nes[tf_name][1]:
                    tf_max_nes[tf_name] = (auc, nes, target_genes )
            
            # Add the max NES values for this run to the global dictionaries
            for tf_name, (auc, nes, target_genes) in tf_max_nes.items():
                tfs_aucs[tf_name].append(auc)
                tfs_nes[tf_name].append(nes)
                tfs_target_genes[tf_name] += (target_genes)

    # Count the number of appearances for each TF
    tfs_appearances = {tf: len(aucs) for tf, aucs in tfs_aucs.items()}
    # Filter TFs that appear in at least 80% of the runs
    consistent_tfs = {tf: count for tf, count in tfs_appearances.items() if count >= 0.01*number_of_repetitions}
    regulons_apperances_of_sign_tfs = defaultdict(list)
    for i in consistent_tfs:
        print(consistent_tfs[i])
        regulons_apperances_of_sign_tfs[i] = Counter(tfs_target_genes[i])
    print(regulons_apperances_of_sign_tfs)
    with open('consistent_tfs_regulons_'+cluster+'_'+condition+'.csv', 'w') as f:
        for key, value in regulons_apperances_of_sign_tfs.items():
            consistent_target_genes = []
            for i in value:
                print(key, tfs_appearances[key])
                if int(value[i]) >= percentage*int(tfs_appearances[key]):
                    consistent_target_genes.append(i)
            
            f.write(str(key)+" : " +str(tfs_appearances[key])+","+ (','.join(consistent_target_genes )) +"\n")

    #print(consistent_tfs)
    consistent_tfs_sorted = (sorted(consistent_tfs.items(), key=lambda x: x[0]) )
    #print(consistent_tfs_sorted)

    consistent_tfs_stats = {}
    for tf in consistent_tfs_sorted:
        mean_auc = float(np.mean([auc for auc in tfs_aucs[tf[0]]]))
        mean_nes = float(np.mean([nes for nes in tfs_nes[tf[0]]]))
        consistent_tfs_stats[tf[0]] = {'appearances': tf[1], 'mean_auc': mean_auc, 'mean_nes': mean_nes}

    # Sort the dictionary by mean_nes (descending)
    sorted_tfs_stats = sorted(
        consistent_tfs_stats.items(),
        key=lambda item: item[1]['mean_nes'],
        reverse=True)

    # Write to CSV
    with open('sorted_tfs_stats'+'_'+cluster+'_'+condition+'.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['TF', 'Appearances', 'Mean AUC', 'Mean NES'])  # header row
        for tf, stats in sorted_tfs_stats:
            writer.writerow([tf, stats['appearances'], stats['mean_auc'], stats['mean_nes']])

    #print("Saved to sorted_tfs_stats.csv")

    #print(consistent_tfs_stats)
    tf_names = [tf for tf, stats in consistent_tfs_sorted]
    mean_nes = [stats['mean_nes'] for tf, stats in consistent_tfs_stats.items()]
    mean_auc = [stats['mean_auc'] for tf, stats in consistent_tfs_stats.items()]

    # Plot histogram
    #plt.figure(figsize=(20, 8))
    #plt.bar(consistent_tfs.keys(), consistent_tfs.values(), color='skyblue')
    #plt.xlabel('Transcription Factors (TFs)', fontsize=12)
    #plt.ylabel('Number of Appearances', fontsize=12)
    #plt.title('Histogram of TFs Appearing More Than 12 Times cluster ', fontsize=14)
    #plt.xticks(rotation=45, ha='right', fontsize=10)
    #plt.tight_layout()
    #plt.show()
    #plt.savefig('tf_histogram'+'_'+cluster+'_'+condition+'.png', dpi=300)  # Save the figure 
    # Create a bar plot for Mean AUC
    #plt.figure(figsize=(20, 8))
    #plt.bar(tf_names, mean_auc, color='lightgreen', label='Mean AUC')
    #plt.xlabel('Transcription Factors (TFs)', fontsize=12)
    #plt.ylabel('Mean AUC', fontsize=12)
    #plt.title('Mean AUC for Consistent TFs', fontsize=14)
    #plt.xticks(rotation=45, ha='right', fontsize=10)
    #plt.tight_layout()
    #plt.show()
    #plt.savefig('mean_auc_histogram'+'_'+cluster+'_'+condition+'.png', dpi=300)  # Save the figure
    # Create a bar plot for Mean NES        
    #plt.figure(figsize=(20, 8))     
    #plt.bar(tf_names, mean_nes, color='lightcoral', label='Mean NES')
    #plt.xlabel('Transcription Factors (TFs)', fontsize=12)
    #plt.ylabel('Mean NES', fontsize=12)     
    #plt.title('Mean NES for Consistent TFs', fontsize=14)
    #plt.xticks(rotation=45, ha='right', fontsize=10)
    #plt.tight_layout()
    #plt.show()
    #plt.savefig('mean_nes_histogram'+'_'+cluster+'_'+condition+'.png', dpi=300)  # Save the figure
    return((consistent_tfs_sorted))


def process_auc_file(auc_filename, number_of_repetitions):
    tfs = []
    for i in range(1, number_of_repetitions+1):
        with open(auc_filename+str(i)+".csv", 'r') as f:
            tfs_list = f.readline()
            tfs += (tfs_list.strip("\n").split(",")[1:])   
    Counter_tfs = Counter(tfs)
    significant_tfs  = ([word for word, occurrences in Counter_tfs.items() if occurrences >= (number_of_repetitions*0.01)])
    print(significant_tfs)

    tf_sums = defaultdict(lambda: 0)
    tf_counts = defaultdict(lambda: 0)
    rownames = None

    for i in range(1, number_of_repetitions+1):
        df = pd.read_csv(f"{auc_filename}{i}.csv", index_col=0)
        df_filtered = df.loc[:, df.columns.intersection(significant_tfs)]
        #print(df_filtered)
        if rownames is None:
            rownames = set(df_filtered.index)
        else:
            rownames &= set(df_filtered.index)

        for tf in df_filtered.columns:
            if tf_counts[tf] == 0:
                tf_sums[tf] = df_filtered.loc[sorted(rownames), tf]
            else:
                tf_sums[tf] += df_filtered.loc[sorted(rownames), tf]
            tf_counts[tf] += 1

    # Step 3: Compute means only where TF was present
    mean_df = pd.DataFrame(index=sorted(rownames))

    for tf in significant_tfs:
        if tf in tf_sums:
            mean_df[tf] = tf_sums[tf] / number_of_repetitions
    mean_df.to_csv(f"{auc_filename}_mean.csv")
    return mean_df


def find_consistent_target_genes(tf):
    target_genes_wt = []
    target_genes_ko = []
    for condition in ["WT", "KO"]:
        for cluster in ['I', 'II', 'III', 'IV', 'V']:
            for regulon_file in ['consistent_tfs_regulons'+'_'+cluster+'_'+condition+'.csv']:
                with open(regulon_file, 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        tf_name = line.split(":")[0]
                        if tf_name == tf and condition=="WT":
                            target_genes_wt += line.split(",")[1::2]
                        if tf_name == tf and condition=="KO":
                            target_genes_ko += line.split(",")[1::2]
    print(f"TF: {tf}, Target Genes KO: {Counter(target_genes_ko)}")
    print(f"TF: {tf}, Target Genes WT: {Counter(target_genes_wt)}")
    print(len(target_genes_wt))
    print(len(target_genes_ko))
    with open('target_gene_occurences_WT_all.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['TF', 'Target Genes'])
        for gene, count in Counter(target_genes_wt).items():
            writer.writerow([tf, gene, count])  
    with open('target_gene_occurences_KO_all.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['TF', 'Target Genes'])
        for gene, count in Counter(target_genes_ko).items():
            writer.writerow([tf, gene, count])
    target_genes_in_both_conditions = (set(target_genes_wt).intersection(target_genes_ko))
    print((target_genes_in_both_conditions))
    with open('common_target_genes_of_'+tf+'_in_both_conditions', 'w') as newfile:
        for gene in target_genes_in_both_conditions:
            newfile.write(gene + "\n") 




def find_consistent_target_genes_between_clusters(tf, cluster):
    target_genes_wt = []
    target_genes_ko = []
    for condition in ["WT", "KO"]:
        for regulon_file in ['consistent_tfs_regulons'+'_'+condition+'_'+cluster+'.csv']:
            with open(regulon_file, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    tf_name = line.split(":")[0]
                    if tf_name == tf and condition=="WT":
                        target_genes_wt += line.split(",")[1::2]
                    if tf_name == tf and condition=="KO":
                        target_genes_ko += line.split(",")[1::2]
    #print(f"TF: {tf}, Target Genes KO: {target_genes_ko}", len(target_genes_ko))
    #print(f"TF: {tf}, Target Genes WT: {target_genes_wt}", len(target_genes_wt))
    #print(f"TF: {tf}, Target Genes WT: {target_genes_wt}", len(target_genes_wt))
    print('KO target genes' , (tf),  len(target_genes_ko))
    print('WT target genes' ,(tf),  len(target_genes_wt))
    with open('target_gene_WT_'+cluster+'_'+tf.strip(' ')+'.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for gene, count in Counter(target_genes_wt).items():
            writer.writerow([gene])  
    with open('target_gene_KO_'+cluster+'_'+tf.strip(' ')+'.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for gene, count in Counter(target_genes_ko).items():
            writer.writerow([gene])
    target_genes_in_both_conditions = (set(target_genes_wt).intersection(target_genes_ko))
    #print((target_genes_in_both_conditions))

    with open('common_target_genes_of_'+tf.strip(' ')+'_in_both_conditions_of_Th1_Th17_'+cluster+'.csv', 'w') as newfile:
        for gene in target_genes_in_both_conditions:
            newfile.write(gene + "\n") 



for condition in ['WT', "KO"]:
    for cluster in ['all']:
        process_regulon_file('regulons_expression_matrix_cluster_Th1_Th17'+condition+'_rep', 50, cluster, condition, percentage=0.3)
        print(process_auc_file('auc_expression_matrix_cluster_Th1_Th17'+condition+'_rep', 50 ))



all_tfs_clean = ['Ube2k', 'Elf1', 'Zfp281', 'Smad3', 'Yy1', 'Ikzf1', 'Bcl11b', 'Bach2', 'Ep300', 'Irf1', 'Dido1', 'Clock', 'Bhlhe40', 'Zbtb1', 'Gabpb1', 'Zfp148', 'Nr2c2', 'Chd1', 'Mbd2', 'Mecp2', 'Zfp407', 'Fli1', 'Gatad1', 'Cnot3', 'Nfkb1', 'E2f2', 'Dnmt1', 'E2f4', 'Rb1', 'Ctcf', 'Tfdp1', 'Elf4', 'Nr3c1', 'Jund', 'Bclaf1', 'Kat7', 'Ahr', 'Tcf12', 'Foxn2', 'Junb', 'Elf2', 'Bach1', 'Ikzf3', 'Sp100']
genes_of_interest_WT = ['Zfp281 ', 'Ahr ', 'Cnot3 ', 'Ube2k ', 'Clock ', 'Ikzf3 ', 'Jund ', 'Bclaf1 ', 'Zbtb1 ', 'Kat7 ', 'Mbd2 ', 'Elf1 ']
genes_of_interest_KO = ['Sp100 ', 'Bhlhe40 ', 'E2f4 ', 'Bach1 ', 'Irf1 ', 'E2f2 ', 'Gabpb1 ', 'Bach2 ', 'Ep300 ', 'Mecp2 ', 'Tcf12 ', 'Gatad1 ', 'Dido1 ', 'Dnmt1 ', 'Elf4 ', 'Chd1 ', 'Zfp407 ' ]


for cluster in ['all']:
    all_tfs_right_annotation = [i+" " for i in all_tfs_clean]
    for tf in all_tfs_right_annotation:
        (find_consistent_target_genes_between_clusters(tf, cluster))




