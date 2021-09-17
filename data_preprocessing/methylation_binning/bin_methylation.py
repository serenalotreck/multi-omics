"""
Calculate bin statistics for methylation data and save as a pickled object.

Author: Serena G. Lotreck 
"""
from os.path import abspath
import argparse 

import pandas as pd 
from collections import defaultdict
import pickle


def get_column_chunks(all_cols, chunksize):
    """
    Yield a generator that has a list of columns to include in a chunk.

    parameters:
        all_cols, list of str: list of column header names 
        chunksize, int: number of columns to read in at a time
    """
    for i in range(0, len(all_cols), chunksize):
        yield all_cols[i:i+chunksize]


def read_methylation_file(filepath):
    """
    Read in a methylation file and transpose for more efficient handling.
    Also splits the locus ID's into two, the seqid and the base pair number.

    Assumes the methylation file is of the shape (n_accessions, n_base_pairs),
    where the column headers are of the form 
        <Chr#>_<bp#>_<methylation_type>_<+/->
    and have no quoted commas.                          

    parameters:
        filepath, str: path to methylation file
                                                                                        
    returns: 
        df, pandas df: methylation df where rows are base pairs and columns
            are accession IDs
    """
    # Read in first row to get the column headers
    print('\nReading file headers...')
    with open(filepath) as f:
        header_str = f.readline()
    headers = header_str.split(',') 

    # Read in chunks and concat to one df
    print('\nReading in file...')
    for i, use_cols in enumerate(get_column_chunks(headers, 100)):
        print(f'Reading chunk {i} of {len(headers)%100}')
        # If it's the first chunk, make this the df 
        if i == 0:
            # Read in 
            df = pd.read_csv(filepath, use_cols=use_cols)
            # Transpose 
            df = df.T
        # Otherwise, transpose and add to the bottom
        else:
            # Read in 
            df_chunk = pd.read_csv(filepath, use_cols=use_cols)
            # Transpose and append  
            df = pd.concat([df, df_chunk.T])

    # Split index into multiindex 
    print('\nSplitting identifiers into chromosome ID and base pair...')
    chr_idx = [df.index.str.split('_')[0]]
    bp_idx = [df.index.str.split('_')[1]]
    df.index = pd.MultiIndex.from_arrays([chr_idx, bp_idx])
    df.rename_axis(("seqid", "bp"))

    print(f'\nDone! Snapshot of dataframe: df.head()')
    return df 


def get_start_row(gff):
    """
    Determine at what row the data starts in a gff file.

    parameters:
        gff, str: path to gff file with genomic features
    
    returns:
        start_row, int: the index of the row where the data starts, 
            which is equivalent to the number of rows to skip
    """
    start_row = 0
    with open(gff) as f:
        line = f.readline()
        line_num = 0 # Keep track of what line we're on 
        while start_row == 0: # Only read lines until we find the start row 
            if line[0] != '#': # If it doesn't start with a comment, it's the first data row 
                start_row = line_num
            else: line_num += 1

    return start_row 


def get_genes(gff):
    """
    Create Gene class instances for genomic features with the type "gene".

    parameters:
        gff, str: path to gff file with genomic features.

    returns:
        genes, dict of list: keys are seqid's (e.g. 'Chr1'), values are
            lists of Gene instances with that seqid
    """
    # Determine at what row the data starts
    start_row = get_start_row(gff)

    # Read in the data 
    gff_df = pd.read_csv(gff, sep='\t', skiprows=start_row,
            names=['seqid', 'source', 'type', 'start', 'end', 'score', 
                    'strand', 'phase', 'attributes'])

    # Filter df to keep only 'gene' rows
    genes_df = gff_df[gff_df['type'] == 'gene']

    # Make gene objects 
    genes = defaultdict(list)
    for row in genes_df.itterrows():
        gene = Gene(row['seqid'], row['start'], row['end'])
        genes[row['seqid']].append(gene)

    return genes 


def main(gff, CG_pres_abs, CHG_pres_abs, CHH_pres_abs, CG_prop, CHG_prop,
        CHH_prop, out_loc, file_prefix):

    # Make gene objects 
    genes = get_genes(gff)

    # Read in methylation datasets 
    methylation_datasets = {'CG_pres_abs' : CG_pres_abs, 
                            'CHG_pres_abs': CHG_pres_abs, 
                            'CHH_pres_abs': CHH_pres_abs, 
                            'CG_prop'     : CG_prop, 
                            'CHG_prop'    : CHG_prop,
                            'CHH_prop'    : CHH_prop}
    methylation_dfs = {}
    for dataset_name, dataset in methylation_datasets.items():
        print(f'\nReading in dataset {dataset_name}')
        methylation_dfs[dataset_name] = read_methylation_file(dataset)

    # Add methylation sites to Gene objects 
    for seqid in genes.keys():
        # Filter methylation datasets
        current_seqid_dfs = {}
        for dataset_name, df in methylation_dfs.items():
            filtered_df = df.loc[dataset_name, :]
            current_seqid_dfs[dataset_name] = filtered_df 
        # Pass to Gene class methods
        genes_to_check = genes[seqid]
        for gene in genes_to_check:
            for dataset_name, seqid_df in current_seqid_dfs.items():
                gene.set_methylation_attr(seqid_dff, dataset_name)

    # Do binning 
    for gene in genes:
        for methylation_type in methylation_datasets.keys():
            gene.set_bin_stat(methylation_type)

    # Save as pickle objects 
    pickle_path = f'{out_loc}/{file_prefix}_binned_methylation_gene_objs.pickle'
    with open(pickle_path, 'wb') as f:
        pickle.dump(genes, f)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Bin methylated sites')

    parser.add_argument('-gff', type=str,
            help='Path to gff3-formatted file with genomic features.')
    parser.add-argument('-CG_pres_abs', type=str,
            help='Path to file with presence/absence data for CG methylation.')
    parser.add-argument('-CHG_pres_abs', type=str,
            help='Path to file with presence/absence data for CHG methylation.')
    parser.add-argument('-CHH_pres_abs', type=str,
            help='Path to file with presence/absence data for CHH methylation.')
    parser.add-argument('-CG_prop', type=str,
            help='Path to file with proportion data for CG methylation.')
    parser.add-argument('-CHG_prop', type=str,
            help='Path to file with proportion data for CHG methylation.')
    parser.add-argument('-CHH_prop', type=str,
            help='Path to file with proportion data for CHH methylation.')
    parser.add_argument('-out_loc', type=str,
            help='Path to directory to save the output.')
    parser.add_argument('-file_prefix', type=str,
            help='String to prepend to output files.')

    args = parser.parse_args()

    args.gff = abspath(args.gff)
    args.CG_pres_abs = abspath(args.CG_pres_abs)
    args.CHG_pres_abs = abspath(args.CHG_pres_abs)
    args.CHH_pres_abs = abspath(args.CHH_pres_abs)
    args.CG_prop = abspath(args.CG_prop)
    args.CHG_prop = abspath(args.CHG_prop)
    args.CHH_prop = abspath(args.CHH_prop)
    args.out_loc = abspath(args.out_loc)

    main(**vars(args))
