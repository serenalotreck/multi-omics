# iRF
## Install the irf
	pip install irf

#### if scikit-learn is needed, make sure that the version of scikit-learn is lower than 0.24, since the presort parameter has been deprecated and removed from the 0.24 scikit-learn. Try the below one:
	pip install scikit-learn==0.23.2

#### For the hpcc user, you can use this: 
	source /mnt/home/peipeiw/Documents/Ath_GS/IRF/IRF/bin/activate
 
## Run the iRF to get the weights of the features for each iteration, and the stability scores of feature combinations
	python IRF_regression.py -df DF_test.csv -df2 Phenotype_value_383_common_accessions_2017_Grimm.csv -test Test.txt -save_name RL_test -y_name RL


