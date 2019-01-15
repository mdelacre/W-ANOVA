# W-ANOVA

1) All simulations were performed using the script "generate nSim random datasets, performing 3 tests on each dataset, extracting p-values and storing them in a file.R"
For 7 distributions, 80 conditions were tested (varying as a function of sd, sd-ratio, sample sizes and sample sizes ratio). 
For each condition, 1.000.000 random datasets were generated, three test were performed on each datasets, and p-values were extracted and stored in .rds files

2) Observed and expected power were computed for each test (F, F* and W) and for each condition, using the script "Observed and expected power for generated sample datasets.R"
==> See "Power.xlsx" for raw observed and expected power.

3) Thanks to te function "Analyses on power.R", tThe three tests (F, F* and W) were compared in terms of:
- raw power
- consistency between observed and expected power ([O-E]/E)  

4) The Type I error rate were computed for each test (F,F* and W) and for each condition, using the script "Type 1 error rate of generated random datasets.R"
==> See "Type I error rate.xlsx"



