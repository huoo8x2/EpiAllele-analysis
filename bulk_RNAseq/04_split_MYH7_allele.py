import pysam
import pandas as pd
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description="Description of your script")
parser.add_argument('--samplename', type=str, required=True, help='sample name')
parser.add_argument('--totalfile', type=str, required=True, help='total file path')
parser.add_argument('--bamfilePath', type=str, required=True, help='bamfile path')
parser.add_argument('--statPath', type=str, required=True, help='stat store path')
parser.add_argument('--snpfile', type=str, required=True, help='snpfile path')

args = parser.parse_args()
outpath = args.bamfilePath
statpath = args.statPath
samplename = args.samplename

# snp info
snp_mp = pd.read_csv(args.snpfile)

# read samfile
samfile_path = args.totalfile
samfile = pysam.AlignmentFile(samfile_path, "rb")

filtered_reads = []
final_groups_list = []
snp_mp[['position_0_base']] = snp_mp[['Location (bp)']]-1
snp_positions = set(snp_mp['position_0_base'])  # Extract SNP positions from the SNP dataframe

for read in samfile:
    read_group = {'WT': 0, 'MUTANT': 0, 'unknown': 0}

    for query_idx, ref_pos in read.get_aligned_pairs():
        if read.mapping_quality == 0:  # filter reads with mapping quality of 0
            continue

        if query_idx is None or ref_pos is None: 
            continue
        
        if ref_pos in snp_positions:
            snp_base_in_query = read.query_sequence[query_idx]  # Get the base in the query sequence
            
            # Classify the read based on the SNP base
            allele = snp_mp.loc[snp_mp['position_0_base'] == ref_pos, 'Observed'].values[0]
            wt_base = snp_mp.loc[snp_mp['position_0_base'] == ref_pos, 'WT'].values[0].split("/")
            mutant_base = snp_mp.loc[snp_mp['position_0_base'] == ref_pos, 'MUTANT'].values[0].split("/")
            
            if snp_base_in_query in wt_base:
                group = 'WT'
            elif snp_base_in_query in mutant_base:
                group = 'MUTANT'
            else:
                group = 'unknown'

            # Update the group info
            read_group[group] += 1

            select_read = read

            # Save SNP-related information
            filtered_reads.append({
                'read_id': read.query_name,
                'ref_read_start': read.reference_start,
                'ref_read_end': read.reference_end,
                'snp_pos': ref_pos+1,
                'snp_allele': allele,
                'WT_base': "/".join(wt_base),
                'MUTANT_base': "/".join(mutant_base),
                'snp_base': snp_base_in_query,
                'group': group
            })
        
    # Determine the final group based on the counts of parental SNP positions
    if read_group['WT'] == 0 and read_group['MUTANT'] == 0 and read_group['unknown'] == 0:
         continue

    if read_group['WT'] == 0 and read_group['MUTANT'] == 0:
        final_group = 'unknown'
    elif read_group['WT'] > read_group['MUTANT']:
        final_group = 'WT'
    elif read_group['WT'] < read_group['MUTANT']:
         final_group = "MUTANT"
    else:
         final_group = 'unknown'

    # save the final parental grouping result
    final_groups_list.append({
                'read_id': read.query_name,
                'reverse_or_not': read.is_reverse,
                'group': final_group
            })


samfile.close()  # Close the samfile

# Save the result to a DataFrame
filtered_reads_df = pd.DataFrame(filtered_reads)
print(filtered_reads_df.group.value_counts())
filtered_reads_df.to_csv(f"{statpath}/{samplename}.read.snp.info.csv", index=False)


final_groups_df = pd.DataFrame(final_groups_list)
print(final_groups_df.group.value_counts())
final_groups_df.to_csv(f"{statpath}/{samplename}.filtered.reads.group.info.csv", index=False)