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
#snp_mp = pd.read_csv("/data02/hukaijie/EpiAllele/data/myh6_snp/mouse_phenome_SNP.csv")
snp_mp = pd.read_csv(args.snpfile)

# read samfile
samfile_path = args.totalfile
samfile = pysam.AlignmentFile(samfile_path, "rb")

filtered_reads = []
final_groups_list = []
snp_mp[['position_0_base']] = snp_mp[['Location (bp)']]-1
snp_positions = set(snp_mp['position_0_base'])  # Extract SNP positions from the SNP dataframe

dba_reads = pysam.AlignmentFile(f"{outpath}/{samplename}.dba.bam", "wb", template=samfile)
c57_reads = pysam.AlignmentFile(f"{outpath}/{samplename}.c57.bam", "wb", template=samfile)

# Filter and classify reads containing SNP positions
for read in samfile:
    read_group = {'DBA': 0, 'C57': 0, 'unknown': 0}

    for query_idx, ref_pos in read.get_aligned_pairs():
        if read.mapping_quality == 0:  # filter reads with mapping quality of 0
            continue

        if query_idx is None or ref_pos is None: 
            continue
        
        if ref_pos in snp_positions:  # Check if the reference position is a SNP
            snp_base_in_query = read.query_sequence[query_idx]  # Get the base in the query sequence
            
            # if reverse complement is needed
            # if read.is_reverse:
            #     reverse = 1
            #     snp_base_in_query = reverse_complement(snp_base_in_query)
            # else:
            #     reverse = 0
            
            # Classify the read based on the SNP base
            allele = snp_mp.loc[snp_mp['position_0_base'] == ref_pos, 'Observed'].values[0]
            dba_base = snp_mp.loc[snp_mp['position_0_base'] == ref_pos, 'DBA/2J'].values[0].split("/")
            c57_base = snp_mp.loc[snp_mp['position_0_base'] == ref_pos, 'C57BL/6J'].values[0].split("/")
            
            if snp_base_in_query in dba_base:
                group = 'DBA'
            elif snp_base_in_query in c57_base:
                group = 'C57'
            else:
                group = 'unknown'

            # Update the group info
            read_group[group] += 1

            # Save SNP-related information
            filtered_reads.append({
                'read_id': read.query_name,
                'ref_read_start': read.reference_start,
                'ref_read_end': read.reference_end,
                'snp_pos': ref_pos+1,
                'snp_allele': allele,
                'DBA_base': "/".join(dba_base),
                'C57_base': "/".join(c57_base),
                'snp_base': snp_base_in_query,
                'group': group
            })
        
    # Determine the final group based on the counts of parental SNP positions
    if read_group['DBA'] == 0 and read_group['C57'] == 0 and read_group['unknown'] == 0:
         continue

    if read_group['DBA'] == 0 and read_group['C57'] == 0:
        final_group = 'unknown'
    elif read_group['DBA'] > read_group['C57']:
        final_group = 'DBA'
    elif read_group['DBA'] < read_group['C57']:
         final_group = "C57"
    else:
         final_group = 'unknown'

    # save the final parental grouping result
    final_groups_list.append({
                'read_id': read.query_name,
                'reverse_or_not': read.is_reverse,
                'group': final_group
            })

    # Write the read to the corresponding SAM file based on parental classification
    if final_group == 'C57':
            c57_reads.write(read)
    elif final_group == 'DBA':
            dba_reads.write(read)

samfile.close()  # Close the samfile
c57_reads.close()
dba_reads.close()

# Save the result to a DataFrame
filtered_reads_df = pd.DataFrame(filtered_reads)
print(filtered_reads_df.group.value_counts())
filtered_reads_df.to_csv(f"{statpath}/{samplename}.read.snp.info.csv", index=False)


final_groups_df = pd.DataFrame(final_groups_list)
print(final_groups_df.group.value_counts())
final_groups_df.to_csv(f"{statpath}/{samplename}.filtered.reads.group.info.csv", index=False)


