
# ABBA 

The R code implements the Augmented Binary Method for Basket Trials.
Files contained in this repository can be used to reproduce the numerical results reported in the paper entitled
Cherlin S and Wason JMS (2024). Augmented Binary Method for Basket Trials (ABBA).

## Basket trial analysis

simulate.R simulates a data set and stores it into a "datalist" object (saved as "datalist.RData)
abba.R analyses the "datalist" object using the ABBA model and stores the resutls in the "res" object (saved as "resAbba.RData")
bin.R analyses the "datalist" object using the BIN model and stores the resutls in the "res" object (saved as "resBin.RData")
prob_success_abba.R computes the log odds ratios based on the probability of success for the ABBA model.
prob_success_bin.R computes the log odds ratios based on the probability of success for the BIN model.

## Stratified analysis

split.R splits the data set into subtrials. Output: datalist1.RData, datalist2.RData, datalist3.RData
abba_strat.R: analyses the "datalist1" object (subtrial1) using the ABBA model and stores the resutls in the "res" object (saved as "resAbba1.RData")
bin_strat.R: analyses the "datalist1" object (subtrial1) using the BINmodel and stores the resutls in the "res" object (saved as "resBin1.RData")
prob_success_abba_strat.R computes the log odds ratios based on the probability of success for the ABBA model for stratified analysis.
prob_success_bin_strat.R computes the log odds ratios based on the probability of success for the BIN model for stratified analysis.

