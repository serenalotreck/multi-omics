"""
Calculate bin statistics for methylation data and save as a pickled object.

Author: Serena G. Lotreck 
"""
from os.path import abspath
import argparse 

from gene import Gene

import pandas as pd 
import datatable as dt
import pickle


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
            else: start_row += 1

    return start_row 


def get_genes(gff):
    """
    Create Gene class instances for genomic features with the type "gene".

    parameters:
        gff, str: path to gff file with genomic features.

    returns:
        genes, list: list of gene objects
    """
    # Determine at what row the data starts
    print('Getting start row...')
    start_row = get_start_row(gff)

    # Read in the data 
    print('Reading in genome data...')
    gff_df = pd.read_csv(gff, sep='\t', skiprows=start_row,
            names=['seqid', 'source', 'type', 'start', 'end', 'score', 
                    'strand', 'phase', 'attributes'])

    # Filter df to keep only 'gene' rows
    print('Getting genes...')
    genes_df = gff_df[gff_df['type'] == 'gene']
    
    # Make gene objects 
    print('Making gene objects...')
    genes = []
    for index, row in genes_df.iterrows():
        gene_name = row.at['attributes'][row.at['attributes'].index('Name=')+5:]
        gene = Gene(row.at['seqid'], row.at['start'], row.at['end'], 
                row.at['strand'], gene_name)
        genes.append(gene)

    return genes 


def main(gff, gbb, ppb, CG_pres_abs, CHG_pres_abs, CHH_pres_abs, CG_prop, 
        CHG_prop, CHH_prop, out_loc, file_prefix):

    # Make gene objects 
    print('\nRetrieving  genes from genome...')
    genes = get_genes(gff)

    # Read in methylation datasets 
    print('\nReading in methylation datasets...')
    methylation_datasets = {'CG_pres_abs' : CG_pres_abs, 
                            'CHG_pres_abs': CHG_pres_abs, 
                            'CHH_pres_abs': CHH_pres_abs, 
                            'CG_prop'     : CG_prop, 
                            'CHG_prop'    : CHG_prop,
                            'CHH_prop'    : CHH_prop}
    
    for dset_name, dset in methylation_datasets.items():
        print(f'Reading in {dset_name}...') 
        # Read in dataset with datatable and make into a pandas df 
        dset_df = dt.fread(dset).to_pandas().set_index('C0').T ## TODO make sure this works 
        print(dset_df.head())
        print('Making multiindex...')
        chr_idx = [idx[0] for idx in dset_df.index.str.split('_')]
        bp_idx = [int(idx[1]) for idx in dset_df.index.str.split('_')]
        strand_idx = [idx[3] for idx in dset_df.index.str.split('_')]
        dset_df.index = pd.MultiIndex.from_arrays([chr_idx, bp_idx, strand_idx],
                                            names=('seqid','bp','strand'))
        # Assign back to dict 
        methylation_datasets[dset_name] = dset_df
    
    # Assign methylation attributes and do binning 
    print('\nBinning data...')
    for i, gene in enumerate(genes):
        print(f'Working on gene {i} of {len(genes)}')
        for dset_name, dset in methylation_datasets.items():
           gene.set_methylation_attr(dset, dset_name)
           gene.set_bin_stat(dset_name, gbb, ppb)

    # Save as pickle objects 
    pickle_path = f'{out_loc}/{file_prefix}_binned_methylation_gene_objs.pickle'
    with open(pickle_path, 'wb') as f:
        pickle.dump(genes, f)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Bin methylated sites')

    parser.add_argument('-gff', type=str,
            help='Path to gff3-formatted file with genomic features.')
    parser.add_argument('-gbb', type=int, default=20,
            help='Number of bins to use for the gene body. Default is 20.')
    parser.add_argument('-ppb', type=int, default=5,
            help='Number of bins to use for pre and post seqs. Default is 5.')
    parser.add_argument('-CG_pres_abs', type=str,
            help='Path to file with presence/absence data for CG methylation.')
    parser.add_argument('-CHG_pres_abs', type=str,
            help='Path to file with presence/absence data for CHG methylation.')
    parser.add_argument('-CHH_pres_abs', type=str,
            help='Path to file with presence/absence data for CHH methylation.')
    parser.add_argument('-CG_prop', type=str,
            help='Path to file with proportion data for CG methylation.')
    parser.add_argument('-CHG_prop', type=str,
            help='Path to file with proportion data for CHG methylation.')
    parser.add_argument('-CHH_prop', type=str,
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
