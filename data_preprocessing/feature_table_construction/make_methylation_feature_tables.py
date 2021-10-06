"""
Use Gene objects containing methylation data to build methylation feature tables.

Author: Serena G. Lotreck 
"""
import argparse
from os.path import abspath
import sys 

sys.path.append('../methylation_binning/')
import gene
import pandas as pd 
import pickle


def make_new_cols(df):
    """
    Make mapper for renaming bins with abbreviated names.
    
    parameters:
        df, pandas df: any one of the 6 methylation stat attribute df's

    returns:
        new_cols, dict: keys are old col names and values are corresponding new 
            col names
    """
    new_cols = {}
    
    for col in df.columns.values.tolist():
        
        bin_num = col.split('_')[-1]
        
        if 'pre' in col:
            new_col = f'R_{bin_num}'
        elif 'gene' in col:
            new_col = f'B_{bin_num}'
        elif 'post' in col:
            new_col = f'P_{bin_num}'
        
        new_cols[col] = new_col

    return new_cols


def make_feature_matrices(genes):
    """
    Make feature matrices for all 6 datasets.

    parameters:
        genes, generator of Gene obj: genes to use to make feature matrices 

    returns:
        data, dict: keys are attribute names, values are full feature matrices 
    """
    # Initialize variables for feature tables
    print('Initializing loop variables...')
    data = {'CG_pres_abs_mean' :None,
            'CG_prop_median'   :None,
            'CHG_pres_abs_mean':None,
            'CHG_prop_median'  :None,
            'CHH_pres_abs_mean':None,
            'CHH_prop_median'  :None}

    # Loop over genes and add to tables 
    print('Looping over genes to build feature matrices...')
    for gene in genes:
        print('my name is elder brown')

        # Check and reverse strands if necessary
        gene.reverse_datasets()

        # Make a dict for renaming columns 
        new_cols = make_new_cols(gene.CG_pres_abs_mean) 
        
        # If it's the first one, make it the df 
        if i == 0:

            for attr_name, value in data.items():
                # Get the df 
                df = getattr(gene, attr_name)
                # Rename cols with shorthand names 
                df = df.rename(columns=new_cols)
                # Add prefix
                df = df.add_prefix(f'{gene.name}_{attr_name}_')
                # Assign back to dict 
                data[attr_name] = df 

        # Otherwise, append
        else:

            for attr_name, value in data.items():
                # Get the new addition df 
                df = getattr(gene, attr_name)
                # Rename cols with shorthand names 
                df = df.rename(columns=new_cols)
                # Add prefix
                df = df.add_prefix(f'{gene.name}_{attr_name}_')
                # Concat onto existing df 
                df = pd.concat(data[attr_name], df, axis=1)
                # Put back in dict 
                data[attr_name] = df

    return data 


def loadall(filename):
    """
    Loads pickled objects as a generator.

    Taken from https://stackoverflow.com/a/28745948/13340814
    """
    with open(filename, "rb") as f:
        while True:
            try:
                print('getting that  pickle')
                yield pickle.load(f)
            except EOFError:
                print('exception reached')
                break

            
def main(pickles, out_loc, file_prefix):

    # Load pickles 
    print('\nLoading genes...')
    genes = loadall(pickles)

    # Make feature matrices 
    print('\nMaking feature matrices...')
    data = make_feature_matrices(genes)

    # Write out as files 
    print('\nWriting out files...')
    for name, matrix in data.items():
        matrix.to_csv(f'{file_prefix}_{name}_feature_matrix.csv') 

    print('\nDone!')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description='Make methylation feature tables.')

    parser.add_argument('pickles', type=str, 
            help='Path to file with Gene object pickles')
    parser.add_argument('out_loc', type=str,
            help='Path to directory to save files')
    parser.add_argument('file_prefix', type=str,
            help='Name to prepend to saved files')

    args = parser.parse_args()

    args.pickles = abspath(args.pickles)
    args.out_loc = abspath(args.out_loc)

    main(args.pickles, args.out_loc, args.file_prefix)
