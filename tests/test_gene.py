"""
Spot checks for the Gene class. 

Author: Serena G. Lotreck 
"""
import unittest 
import sys
import os
import shutil

sys.path.append('../data_preprocessing/methylation_binning/')
from gene import Gene

import pandas as pd
from pandas.testing import assert_frame_equal


class TestGene(unittest.TestCase):

    def setUp(self):

        # Create test directory 
        self.tmpdir = "tmp"
        os.makedirs(self.tmpdir, exist_ok=False)
        
        # Create test files 
        CG_pres_abs = (',Chr1_23_CG_+,Chr1_553_CG_+,Chr1_872_CG_-,Chr1_1230_CG_+,'
                'Chr2_720_CG_+\nac1,1,0,1,1,0\nac2,'
                '0,0,1,0,1')
        self.CG_pres_abs = f'{self.tmpdir}/CG_pres_abs.csv'
        with open(self.CG_pres_abs, 'w') as myf:
            myf.write(CG_pres_abs)

        CG_prop = (',Chr1_23_CG_+,Chr1_553_CG_+,Chr1_872_CG_-,Chr1_1230_CG_+,'
                'Chr2_720_CG_+\nac1,0.23,0.45,0.21,0.97,0.11\nac2,'
                '0.65,0.23,0.67,0.19,0.87')   
        self.CG_prop = f'{self.tmpdir}/CG_prop.csv'
        with open(self.CG_prop, 'w') as myf:
            myf.write(CG_prop)
        
        # Create the Gene instance
        self.gene = Gene('Chr1', 560, 890, '+', 'AT102947')


    def tearDown(self):

        shutil.rmtree(self.tmpdir)

        
    def test_set_methylation_attr(self):

        self.gene.set_methylation_attr(self.CG_pres_abs, 'CG_pres_abs')

        right_answer = pd.DataFrame({'ac1':[0,1,1], 'ac2':[0,1,0]},
                index=pd.MultiIndex.from_arrays([['Chr1','Chr1','Chr1'],
                                                 [553, 872, 1230],
                                                 ['+', '-', '+']],
                                                names=('seqid', 'bp', 'strand')))
        
        assert_frame_equal(self.gene.CG_pres_abs, right_answer)


    def test_bin_data_calculate_stat_pre_mean(self):

        self.gene.set_methylation_attr(self.CG_pres_abs, 'CG_pres_abs')
        
        pre_df = self.gene.CG_pres_abs.loc[(self.gene.start > 
            self.gene.CG_pres_abs.index.get_level_values('bp')) & 
            (self.gene.start-500 < 
            self.gene.CG_pres_abs.index.get_level_values('bp'))]
        
        pre_stat_df_mean = self.gene.bin_data_calculate_stat(pre_df,
                self.gene.start-500, self.gene.start, 'CG_pres_abs', 'pre', 5)
        
        right_answer = pd.DataFrame({f'pre_bin_{i}':[0.0, 0.0] for i in range(5)},
                index=['ac1','ac2'])

        # Don't check dtype because if mean or median functions are used on only 
        # one item, if that item is an int, the result is an int
        assert_frame_equal(pre_stat_df_mean, right_answer, check_dtype=False)


    def test_bin_data_calculate_stat_gene_mean(self):
       
        self.gene.set_methylation_attr(self.CG_pres_abs, 'CG_pres_abs')
        
        gene_body_df = self.gene.CG_pres_abs.loc[(self.gene.start < 
            self.gene.CG_pres_abs.index.get_level_values('bp')) & 
            (self.gene.CG_pres_abs.index.get_level_values('bp') < 
            self.gene.end)]
        
        gene_stat_df_mean = self.gene.bin_data_calculate_stat(gene_body_df,
                self.gene.start, self.gene.end, 'CG_pres_abs', 'gene_body', 20)


        right_ans_dicts = {**{f'gene_body_bin_{i}':[0.0,0.0] for i in range(18)},
                        **{'gene_body_bin_18':[1.0,1.0]},
                        **{'gene_body_bin_19':[0.0,0.0]}}
        right_answer = pd.DataFrame(right_ans_dicts, index=['ac1','ac2'])
            
        assert_frame_equal(gene_stat_df_mean, right_answer, check_dtype=False)
    
    
    def test_bin_data_calculate_stat_post_median(self):

        self.gene.set_methylation_attr(self.CG_prop, 'CG_prop')
        
        post_df = self.gene.CG_prop.loc[(self.gene.end < 
            self.gene.CG_prop.index.get_level_values('bp')) 
            & (self.gene.end+500 > self.gene.CG_prop.index.get_level_values('bp'))]
        
        post_stat_df_median = self.gene.bin_data_calculate_stat(post_df,
                self.gene.end, self.gene.end+500, 'CG_prop', 'post', 5)


        right_ans_dict = {**{f'post_bin_{i}':[0.0,0.0] for i in range(3)},
                          **{'post_bin_3':[0.97, 0.19]},
                          **{'post_bin_4':[0.0,0.0]}}
        right_answer = pd.DataFrame(right_ans_dict, index=['ac1','ac2'])

        assert_frame_equal(post_stat_df_median, right_answer, check_dtype=False)
    
    
    def test_set_bin_stat_mean(self):

        pass

    def test_bin_stat_median(self):
    
        pass

if __name__ == '__main__':
    unittest.main()

                                        
