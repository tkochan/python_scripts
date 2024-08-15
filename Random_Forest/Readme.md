Random Forest: 

Take Clustage output files:
1.	out_subelements.csv
2.	out_subelements.key.txt

1.	Variable Selection:
a.	Make three Models:
i.	Unfiltered
1.	Split strains into training (n= 207) and test (n=50), use all variables
2.	Train a model on all variables and evaluate model fit for the test dataset
3.	random_forest_tuning_unfiltered.py

ii.	Filtered
1.	Split strains into training (n= 207) and test (n=50), and filter variables based on length, num_genomes, and avg_rank.
2.	Train a model on the filtered variables and evaluate model fit for the test dataset
3.	random_forest_tuning.py


iii.	Filtered & Weighted
1.	Split strains into training (n= 207) and test (n=50), and filter variables based on length, num_genomes, and avg_rank.
2.	Split the training further into 182 vs 25 at random 100x.
3.	Iterate through each random split 30x and add variables 1 at a time choosing the variable that has the greatest increase in KSF (Kochan statistical factor).  
4.	Combine all variables into a CSV and sort by the variables that have the highest number of occurrences. 
5.	Train a model using the 50 – 250 top variables in increments of 5 (with original train / test split) and choose the model that has the lowest test MSE.  
6.	Random_forest_weighted_variable_selection.py (fix this to output the variables and number of occurrences)
7.	random_forest_tuning_weighted3.py

iv.	Filtered and Efron-Gong Bootstrap Optimization
1.	Split strains into training (n= 207) and test (n=50), and filter variables based on length, num_genomes, and avg_rank.
2.	Bootstrap samples 10x
3.	Iterate through variables 1 at a time (30x) performing regression model on each bootstrap. Substract average bootstrap_r2 from training_r2 to determine optimism. Choose the variable that has the highest AdjR2 (R2-optimism)
4.	Train a model using the 30 variables.  
5.	
6.	test.py (fix this to output the variables and number of occurrences)
7.	random_forest_tuning_weighted3.py


b.	Hyperparameter tuning
i.	Unfiltered
1.	Split strains into training (n= 207) and test (n=50), use all variables 
2.	Use a grid search and 5x fold cross validation on training data set to evaluate performance of different hyperparameter combinations
3.	random_forest_tuning_unfiltered.py


ii.	Filtered
1.	Split strains into training (n= 207) and test (n=50), and filter variables based on length, num_genomes, and avg_rank.
2.	Use a grid search and 2x fold cross validation on training data set to evaluate performance of different hyperparameter combinations
3.	random_forest_tuning.py


iii.	Filtered & Weighted
1.	Split strains into training (n= 207) and test (n=50), and use top 70 variables.
2.	Use a grid search and 2x fold cross validation on training data set to evaluate performance of different hyperparameter combinations.
3.	random_forest_tuning_weighted.py

iv.	Save Hyperparameter settings
1.	Save hyperparameter settings for each model in the file: hyperparemeters_tuning.csv



c.	Model Comparisons
i.	Unfiltered
1.	Split strains into training (n= 207) and test (n=50), use all variables and best combination of hyperparameters
2.	Train a model on all variables and evaluate model fit for the test dataset

ii.	Filtered
1.	Split strains into training (n= 207) and test (n=50), use filtered variables and best combination of hyperparameters
2.	Train a model on filtered variables and evaluate model fit for the test dataset 
3.	
iii.	Filtered & Weighted
1.	Split strains into training (n= 207) and test (n=50), and use top 70 variables and best combination of hyperparameters
2.	Train a model on top 70 variables and evaluate model fit for the test dataset 
iv.	Compare Models
1.	Compare training r2 and test r2 for each model
2.	Compare OTHER statistics
3.	Determine accuracy and change output to binary hmv vs non-hmv
4.	Create ROC curve and determine AUC
5.	random_forest_comparisons.py









Kochan Statistical Factor (KSF)

ksf = r2_train + 1.2 * r2_test

The idea is that out of thousands of variables, you want to choose the variables that have the most general applicability. If you only choose the variables that improve the r2_train you are choosing the variables that have the greatest correlation with the training data set. This is a recipe for over fitting and poor performance with a future test data set. 

By further splitting the training data set into train and test data sets randomly 100x and then using ksf, the model will select the variables that have the highest r2_train and also the higher r2_test with r2_test being valued with an extra 20%. 

In this way, if the top two variables have an r2_train of 0.8 it will choose the variable that has the higher r2_test. 

This way variables are selected that when trained on a training data set they produce accurate predictions on the test data set. This would theoretically create a more model that produce accurate predictions for a dataset the model has never seen before.  
