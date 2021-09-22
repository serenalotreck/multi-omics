"""
Calculate bin statistics for methylation data and save as a pickled object.

Author: Serena G. Lotreck 
"""
from os.path import abspath
import argparse 

from gene import Gene

import pandas as pd 
import pickle


def check_column_constistency(filepath):
    """
    Check that all rows have the same number of columns when split by 
    commas. 

    parameters:
        filepath, str: path to file to check 

    returns: True if all rows are the same, False otherwise
    """
    with open(filepath) as f:
        lines = readlines()
    lines = [line.split(',') for line in lines]

    if len({len(line) for line in lines}) == 1:
        return True 
    else:
        return False 


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


def main(gff, CG_pres_abs, CHG_pres_abs, CHH_pres_abs, CG_prop, CHG_prop,
        CHH_prop, out_loc, file_prefix):

    # Make gene objects 
    print('\nRetrieving  genes from genome...')
    genes = get_genes(gff)

    # Read in methylation datasets 
    methylation_datasets = {'CG_pres_abs' : CG_pres_abs, 
                            'CHG_pres_abs': CHG_pres_abs, 
                            'CHH_pres_abs': CHH_pres_abs, 
                            'CG_prop'     : CG_prop, 
                            'CHG_prop'    : CHG_prop,
                            'CHH_prop'    : CHH_prop}
    
    # Check that the files are all valid 
    for dset_name, dset in methylation_datasets.items():
        assert check_column_constistency(dset), (
                f'There is at least one row in {dset_name} that has '
                'a different number of columns than the rest of the data. '
                'Please fix before continuing.')
    
    # Assign methylation attributes and do binning 
    for gene in genes:
        for dset_name, dset in methylation_datasets.items():
           gene.set_methylation_attr(dset, dset_name)
           gene.set_bin_stat(dset_name)

    # Save as pickle objects 
    pickle_path = f'{out_loc}/{file_prefix}_binned_methylation_gene_objs.pickle'
    with open(pickle_path, 'wb') as f:
        pickle.dump(genes, f)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Bin methylated sites')

    parser.add_argument('-gff', type=str,
            help='Path to gff3-formatted file with genomic features.')
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
