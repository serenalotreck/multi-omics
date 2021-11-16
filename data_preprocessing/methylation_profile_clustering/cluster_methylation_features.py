"""
Make a feature table based on clustering of methylation profiles.

Author: Serena G. Lotreck
"""
import argparse
from os.path import abspath

import datatable as dt
from sklearn.cluster import KMeans


def main(feature_table, num_clusters, out_loc):

    # Read in data
    print('\nReading in feature table...\n')
    df = dt.fread(feature_table).to_pandas.set_index('C0').T
    print(df.head())

    # Separate into gene-level bins
    ## come back to this


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Cluster methylation data')

    parser.add_argument('feature_table', type=str,
            help='Path to feature table to cluster')
    parser.add_argument('num_clusters', type-int,
            help='Number of clusters to use for kmeans')
    parser.add_argument('out_loc', type=str,
            help='Path to save files')

    args = parser.parse_args()

    args.feature_table = abspath(args.feature_table)
    args.out_loc = abspath(args.out_loc)

    main(args.feature_table, args.num_clusters, args.out_loc)
