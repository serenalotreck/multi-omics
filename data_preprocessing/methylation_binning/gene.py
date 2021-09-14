"""
Gene class definition to use for binning methylation.

Author: Serena G. Lotreck
"""

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

        returns: None 
        """
        # Check datatype of index & change to int if it's not int
        if methylation_df.index.dtype != 'int64':
            methylation_df.index = methylation_df.index.astype('int64')
            
        # Filter df by whether or not it's in the gene 
        in_near_gene_df = methylation_df.loc[(self.start-500 < methylation_df.index) & 
                                            (methylation_df.index < self.end+500)]

        # Set attribute 
        setttr(self, attr_name, in_near_gene_df)
    

    def set_bin_stat(self, methylation_type):
        """
        Separate sites into bins and calculate statistic across bins and 
        accessions for a given type of methylation, and assign the result 
        as an attribute.

        NOTE: If there are no C's in a bin, there will be a NaN, as this 
        seems to more accurately reflect what we're trying to demonstrate
        rather than putting a 0 there. This may become a problem downstream,
        which is why I mention it now. 

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
        upstream_df = df.loc[~(self.start < df.index)]
        downstream_df = df.loc[~(df.index < self.end)]

        # Bin and calculate within gene bins
        
        ## Within gene body 
        min_bp = df.bp.min()
        max_bp = df.bp.max()
        num_gene_bins = 20
        gene_body_df.bp = pd.cut(gene_body_df.bp, bins=np.arange(min_bp, max_bp, 
                    step=(max_bp - min_bp)/num_gene_bins, include_lowest=True))
        if 'pres_abs' in methylation_type:
            gene_stat_df = gene_body_df.groupby('bp').mean()
        elif 'prop' in methylation_type:
            gene_stat_df = gene_body_df.groupby('bp').median()
        # Rename bins ##TODO make this work 
        bins = 
        gene_bins = {bn:f'gene_body_bin{i}' for i, bn in enumerate(bins)}
        gene_stat_df = gene_stat_df.T.rename(columns=gene_bins)


        ## Upstream 
        start_range = self.start - 500
        num_up_bins = 5 # 100bp bins
        upstream_df.bp = pd.cut(upstream_df.bp, bins=np.arange(start_range, 
                    self.start, step=(self.start - start_range)/num_up_bins,
                    include_lowest=True)) 
        if 'pres_abs' in methylation_type:
            up_stat_df = upstream_df.groupby('bp').mean()
        elif 'prop' in methylation_type:
            up_stat_df = upstream_df.groupby('bp').median()
        # Rename bins
        ##TODO
        
        ## Downstream 
        end_range = self.end + 500
        num_down_bins = 5 # 100bp bins
        downstream_df.bp = pd.cut(downstream_df.bp, bins=np.arange(self.end, 
                    end_range, step=(end_range - self.end)/num_down_bins,
                    include_lowest=True)) 
        if 'pres_abs' in methylation_type:
            down_stat_df = downstream_df.groupby('bp').mean()
        elif 'prop' in methylation_type:
            down_stat_df = downstream_df.groupby('bp').median()
        ##TODO
