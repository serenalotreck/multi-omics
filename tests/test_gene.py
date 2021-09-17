"""
Spot checks for the Gene class. 

Author: Serena G. Lotreck 
"""
import unittest 
import sys
sys.path.append('../data_preprocessing/methylation_binning/')
from gene import Gene
import pandas as pd
from pandas.testing import assert_frame_equal


class TestGene(unittest.TestCase):

    def setUp(self):

        self.CG_pres_abs = pd.DataFrame({'accession1':[1,0,1,1], 
                                         'accession2':[0,0,1,0]},
                                        index=[23,553,872,1230])
        self.CG_prop = pd.DataFrame({'accession1':[0.23,0.45,0.21,0.97],
                                     'accesion2':[0.65,0.23,0.67,0.19]},
                                        index=[23,553,872,1230])
        self.gene = Gene('Chr1', 560, 890)

        # Make stat df's to check answers 
        ## Column headers 
        up_bin_names = [f'upstream_bin_{i}' for i in range(5)]
        down_bin_names = [f'downstream_bin_{i}' for i in range(5)]
        gene_bin_names = [f'gene_body_bin_{i}' for i in range(20)]
        columns = up_bin_names + gene_bin_names + down_bin_names
        ## Data
        up_mean = [float(i) for i in [0, 0, 0, 0, 0]]
        down_mean = [float(i) for i in [0, 0, 0, 0.5, 0]]
        gene_mean = [float(i) for i in [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]]
        up_median = [float(i) for i in [0, 0, 0, 0, 0.34]]
        down_median = [float(i) for i in [0, 0, 0, 0.58, 0]]
        gene_median = [float(i) for i in [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.44, 0]]
        ## Dataframes 
        self.setattr_CG_pres_abs = pd.DataFrame({'accession1':[0,1,1],
                                                'accession2':[0,1,0]},
                                                index=[553,872,1230])
        self.up_stat_df_mean = pd.DataFrame([up_mean], columns=up_bin_names)
        self.gene_stat_df_mean = pd.DataFrame([gene_mean], columns=gene_bin_names)
        self.down_stat_df_median = pd.DataFrame([down_median], columns=down_bin_names)
        self.CG_pres_abs_mean = pd.DataFrame([up_mean+gene_mean+down_mean], 
                columns=columns)
        self.CG_prop_median = pd.DataFrame([up_median+gene_median+down_median], 
                columns=columns)


    def test_set_methylation_attr(self):

        self.gene.set_methylation_attr(self.CG_pres_abs, 'CG_pres_abs')
        
        assert_frame_equal(self.gene.CG_pres_abs, self.setattr_CG_pres_abs)


    def test_bin_data_calculate_stat_upstream_mean(self):

        upstream_df = self.setattr_CG_pres_abs.loc[(self.gene.start > 
            self.setattr_CG_pres_abs.index) & (self.gene.start-500 < 
            self.setattr_CG_pres_abs.index)]
        up_stat_df_mean = self.gene.bin_data_calculate_stat(upstream_df,
                self.gene.start-500, self.gene.start, 'CG_pres_abs', 'upstream')
        
        assert_frame_equal(up_stat_df_mean, self.up_stat_df_mean)


    def test_bin_data_calculate_stat_gene_mean(self):
        
        gene_body_df = self.CG_pres_abs.loc[(self.gene.start < 
            self.CG_pres_abs.index) & (self.CG_pres_abs.index < self.gene.end)]
        gene_stat_df_mean = self.gene.bin_data_calculate_stat(gene_body_df,
                self.gene.start, self.gene.end, 'CG_pres_abs', 'gene_body')

        assert_frame_equal(gene_stat_df_mean, self.gene_stat_df_mean)
    
    
    def test_bin_data_calculate_stat_downstream_median(self):

        downstream_df = self.CG_prop.loc[(self.gene.end < 
            self.CG_prop.index) & (self.gene.end+500 > self.CG_prop.index)]
        down_stat_df_median = self.gene.bin_data_calculate_stat(downstream_df,
                self.gene.end, self.gene.end+500, 'CG_prop', 'downstream')

        assert_frame_equal(down_stat_df_median, self.down_stat_df_median)
    
    
    def test_set_bin_stat_mean(self):

        self.gene.set_methylation_attr(self.CG_pres_abs, 'CG_pres_abs')
        self.gene.set_bin_stat("CG_pres_abs")

        assert_frame_equal(self.gene.CG_pres_abs_mean, self.CG_pres_abs_mean)


    def test_bin_stat_median(self):

        self.gene.set_methylation_attr(self.CG_prop, 'CG_prop')
        self.gene.set_bin_stat("CG_prop")

        assert_frame_equal(self.gene.CG_prop_median, self.CG_prop_median)

if __name__ == '__main__':
    unittest.main()

                                        
