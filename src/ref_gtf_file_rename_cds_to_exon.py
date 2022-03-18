#!/usr/bin/env python3
#%%
import pandas as pd
import numpy as np
import csv
import re
import argparse

#%%

def replace_transcript_id_with_pacbio_info(row):
    transcript_id_raw = row['attribute']
    result = re.search(r'transcript_id "(.*\|)(.*)(\|.*)"', transcript_id_raw)
    transcript_id = re.sub(r'transcript_id .*', f'transcript_id "{result.group(2)}";', transcript_id_raw)
    return transcript_id

def process_reference_rename(reference_file, name):
    gtf_colnames = ['seqname','source','feature','start','end','score','strand','frame','attribute']
    ref = pd.read_table(reference_file, names=gtf_colnames, skiprows=5, )
    ref_cds = (
        ref
            .query('feature in ["transcript","CDS"]')
            .copy()  
    )
    ref_cds['feature'] = np.where(ref_cds['feature'] =='CDS','exon','transcript')

    savefile = f'{name}.cds_renamed_exon.gtf'
    # with open(reference_file, 'r') as iref,open(savefile, 'w') as out:
    #     for i in range(5):
    #         out.write(iref.readline())
    ref_cds.to_csv(savefile, sep='\t',index=False, header=False, quoting=csv.QUOTE_NONE)

    ref_exons = (
        ref
            .query('feature in ["transcript","exon"]')
    )
    savefile = f'{name}.transcript_exon_only.gtf'
    # with open(reference_file, 'r') as iref, open(savefile, 'w') as out:
    #     for i in range(5):
    #         out.write(iref.readline())
    ref_exons.to_csv(savefile, sep='\t',index=False, header=False, quoting=csv.QUOTE_NONE)


def main():
    parser = argparse.ArgumentParser(description='rename cds to exon for sqanti protein module')
    parser.add_argument('--reference_gtf', action='store', dest= 'reference_gtf',help='sample gtf file')
    parser.add_argument('--reference_name', action='store', dest= 'reference_name',help='sample name')
    results = parser.parse_args()
    
    process_reference_rename(results.reference_gtf, results.reference_name)

    
if __name__ == "__main__":
    main()

# %%
