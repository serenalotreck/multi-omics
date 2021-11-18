"""
Make a feature table based on clustering of methylation profiles.

Author: Serena G. Lotreck
"""
import argparse
from os.path import abspath, basename
from collections import defaultdict
import itertools

import datatable as dt
from sklearn.cluster import KMeans
from sklearn.preprocessing import OneHotEncoder


def make_tidy_data(df):
    """
    Make a matrix of shape (n_accessions*n_genes, n_bins), where each row
    is {instance (accession) ID}_{gene #}, and the columns are the bins that
    methylation data has been placed into.

    The columns of the input dataset must have a very specific format:
        {gene name}_{methylation & data type}_{R/B/P}_{bin number}
    for example:
        AT1G01010_CG_pres_abs_mean_R_0
    There can be underscores in the methylation and data type, but there
    can NOT be underscores in the gene name, and the bin ID (R/B/P) and number
    must be the last two items when the string is split on underscores.

    parameters:
        df, pandas df: rows are instances and columns are the 30 bins for each
            of 27,000 genes

    returns:
        tidy_df, pandas df: the tidy version of the data
    """
    # Get the lists of column names that correpsond to each gene
    cols = df.columns.values.tolist()
    gene_names = lambda x: x.split('_')[0] # Lambda func for groupby
    col_lists = defaultdict(list)
    for gene_name, col_names in itertools.groupby(cols, key=gene_names):
        col_lists[gene_name] += list(col_names)

    # Make df
    for gene_name, col_names in col_lists.items():

        # Make a mapper for the columns
        bin_names = lambda x: '_'.join(x.split('_')[-2:])
        column_mapper = {name:bin_names(name) for name in col_names}

        # Check if this is the first one
        try:
            tidy_df
        except NameError:
            tidy_df = df[col_names]
            tidy_df.rename(columns=column_mapper)

       # Otherwise just concat them
       ## TODO come back to this


def main(feature_table_path, num_clusters, out_loc):

    # Read in data
    print('\nReading in feature table...\n')
    df = dt.fread(feature_table_path).to_pandas().set_index('C0')
    print(df.head())

    # Preprocess data
    accessions = df.index
    X = df.to_numpy()

    # Cluster
    print('\nClustering...')
    clusters = KMeans(n_clusters=num_clusters, random_state=42).fit_predict(X)

    # One hot encode
    print('\nOne hot encoding...')
    encoder = OneHotEncoder()
    encoded = encoder.fit_transform([clusters])

    # Make feature table
    print('\nMaking feature table...')
    feats = encoder.get_feature_names_out(input_features=['cluster'])
    print(f'Feature names: {feats[:5] if len(feats) > 5 else feats}')
    feature_table = pd.DataFrame(encoded, columns=feats, index=accessions)
    print(f'Snapshot of feature table:\n{feature_table.head()}')

    # Write out
    print('\nWriting out file...')
    file_name = basename(feature_table_path)
    feature_table.to_csv(f'{out_loc}/{file_name}_{num_clusters}_clusters.csv')

    print('\nDone!')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Cluster methylation data')

    parser.add_argument('feature_table', type=str,
            help='Path to feature table to cluster')
    parser.add_argument('num_clusters', type=int,
            help='Number of clusters to use for kmeans')
    parser.add_argument('out_loc', type=str,
            help='Path to save files')

    args = parser.parse_args()

    args.feature_table = abspath(args.feature_table)
    args.out_loc = abspath(args.out_loc)

    main(args.feature_table, args.num_clusters, args.out_loc)
