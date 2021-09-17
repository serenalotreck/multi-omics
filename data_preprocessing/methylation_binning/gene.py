"""
Gene class definition to use for binning methylation.

Author: Serena G. Lotreck
"""
import pandas as pd 
import numpy as np
import operator 
import re


class Gene:

    def __init__(self, seqid, start, end):
        """
        Instantiate an instance of the Gene class. 

        parameters:
            seqid, str: ID from the first column of the gff file 
            start, int: start bp
            end, int: end bp

        returns: None
        """
        self.seqid = seqid
        self.start = start
        self.end = end

        # These will be df's later 
        self.CG_pres_abs = None
        self.CHG_pres_abs = None
        self.CHH_pres_abs = None
        self.CG_prop = None
        self.CHG_prop = None
        self.CHH_prop = None
        self.CG_pres_abs_mean = None
        self.CHG_pres_abs_mean = None
        self.CHH_pres_abs_mean = None
        self.CG_prop_median = None
        self.CHG_prop_median = None
        self.CHH_prop_median = None

        
    def set_methylation_attr(self, methylation_df, attr_name):
        """
        Filter methylation dataset by what sites are in this gene, 
        or within 500bp up- or downstream, and add
        the resulting df to the attribute corresponding to attr_name.

        parameters:
            methylation_df, df: index is the base pair at which the methylation
                occurrs, column headers are accessions
            attr_name, str: name of attribute to set. 

        returns: None 
        """
        # Check datatype of index & change to int if it's not int
        if methylation_df.index.dtype != 'int64':
            methylation_df.index = methylation_df.index.astype('int64')
            
        # Filter df by whether or not it's in the gene 
        in_near_gene_df = methylation_df.loc[(self.start-500 < methylation_df.index) & 
                                            (methylation_df.index < self.end+500)]

        # Set attribute 
        setattr(self, attr_name, in_near_gene_df)
    

    @staticmethod
    def bin_data_calculate_stat(df, start, end, methylation_type, name_part):
        """
        Put data into base pair bins for a section of the gene (upstream, 
        gene body, downstream). 

        parameters:
            df, pandas df: the data to bin, accessions are columns and bp
                is the index 
            start, int: what base pair is the start of the region to bin
            end, int: what abse pair is the end of the region to bin 
            methylation_type, str: name of methylation type this is for. Name 
                should correspond to one of the methylation dataset attributes.
            name_part, str: "upstream", "gene_body" or "downstream. Determines 
                the number of bins, 5 for up and downstream, 20 for gene body.

        returns: 
            stat_df, df: accessions are rowns, bins are columns
        """
        # Determine bin number 
        num_bins = 20 if name_part == "gene_body" else 5

        # Get the bin boundaries using qcut on all bases in gene 
        bin_intervals, bins = pd.qcut(pd.Series(np.arange(start, end, 1)), 
                                    num_bins, retbins=True)

        # Then use bin boundaries with cut on the methylated bases 
        df.index = pd.cut(df.index, bins=bins, include_lowest=True)
        
        # Add columns of zeros for bins with no C's in them 
        df = df.T
        num_accs = df.shape[0]
        df = pd.DataFrame({col_name:(np.zeros(num_accs, dtype=int) \
                        if col_name not in df.columns.values.tolist() \
                        else df[col_name].tolist()) for col_name in set(bin_intervals)})
       
        # Groupby to get stats
        if 'pres_abs' in methylation_type:
            stat_df = df.mean(axis=0).to_frame().T
        elif 'prop' in methylation_type:
            stat_df = df.median(axis=0).to_frame().T
        
        # Rename bins
        bin_names = sorted(set(bin_intervals), key=operator.attrgetter('left'))
        col_bins = {bn:f'{name_part}_bin_{i}' for i, bn in enumerate(bin_names)}
        stat_df = stat_df.rename(columns=col_bins)
        
        # Sort columns 
        natsort = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)',s)]
        sorted_cols = sorted(stat_df.columns.values.tolist(), key=natsort)
        stat_df = stat_df.reindex(sorted_cols, axis=1)
        
        return stat_df


    def set_bin_stat(self, methylation_type):
        """
        Separate sites into bins and calculate statistic across bins and 
        accessions for a given type of methylation, and assign the result 
        as an attribute.

        NOTE: If there are no C's in a bin, the values will be set to 0. 

        parameters:
            methylation_type, str: name of methylation type to do this for.
                Name should correspond to one of the methylation dataset 
                attributes. 

        returns: None 
        """ 
        # Get the dataset
        df = getattr(self, methylation_type)

        # Separate into in-gene and before/after-gene sites 
        gene_body_df = df.loc[(self.start < df.index) & (df.index < self.end)]
        upstream_df = df.loc[(self.start > df.index) & (self.start-500 < df.index)]
        downstream_df = df.loc[(df.index > self.end) & (self.end+500 > df.index)]

        # Bin and calculate within gene bins
        gene_body_stat_df = Gene.bin_data_calculate_stat(gene_body_df, self.start, 
                self.end, methylation_type, "gene_body") 
        up_stat_df = Gene.bin_data_calculate_stat(upstream_df, self.start-500, 
                self.start, methylation_type, "upstream") 
        down_stat_df = Gene.bin_data_calculate_stat(downstream_df, self.end, 
                self.end+500, methylation_type, "downstream") 

        # Combine all df's 
        overall_stat_df = pd.concat([up_stat_df, gene_body_stat_df, down_stat_df], axis=1)

        # Set as attribute 
        if 'pres_abs' in methylation_type:
            setattr(self, f'{methylation_type}_mean', overall_stat_df)
        elif 'prop' in methylation_type:
            setattr(self, f'{methylation_type}_median', overall_stat_df)
