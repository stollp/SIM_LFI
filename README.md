---
output:
  pdf_document: default
  html_document: default
---
# nncrexm
R-Scripts

# Crown extension and plasticity

This repository contains ...

## File listing

</p> </center>

### 1.) R-Scripts

<center> <p style="text-align:center;">

| File                     | task                                                 |sourcing/calling   |
| ------------------------ | -----------------------------------------------------|-------------------|
| 0_Impressionistic.R      | read data, calculate crown radius, distribute points |BM.distribute.R 0_Prune.R 0_Prune_CALC.R |
0_Prune.R                  | evaluates neighbourhoods of points within DELTA_D    |                   |
0_Prune_CALC.R             | calculates relocations                               |                   |
| 1_Loop.R                 | loop over species to                                 | 2_NN.R            |
| 2_NN.R                   | evaluate neighbourhood of trees (r = 1, 2, ..., 20 m)| eval_nn           |

</p> </left>

### 2.) R-Functions

<center> <p style="text-align:center;">

| File                     | task                                                                        | function name |
| ------------------------ | --------------------------------------------------------------------------- |-------------- |
| BM_Distribute.R          | distribute points                                                           |BM.distribute  |
|       eval_NN.R          | evaluate neigbhourhoods of trees                                            |eval_nn        |

</p> </center>


## Variable listing

<center>

<p style="text-align:center;">

| variable | meaning                      |
| -------- | ---------------------------  |
| Plot     | Main plot 1 (MP1) or 2 (MP2) |
| Tag      | individual tag number        |
|   X      | X-coord in m                 |
|   Y      | Y-coord in m                 |
| F_Code07 | Family code                  |
|   Code07 | Species code                 |
| GBH86    | Girth at breast height 1986  |
| GBH96    | Girth at breast height 1996  |
| GBH01    | Girth at breast height 2001  |
| GBH07    | Girth at breast height 2007  |

</p>

</center>

Data can be downloaded from a [Dryad repository] ...

### Summary of variables in data file

```{r eval=TRUE, echo=TRUE, fig.height=4}
d <- read.table("FULL_test.txt", header = TRUE)
summary(d[, c("Plot", "X", "Y", "Code07", "GBH86", "GBH96", "GBH01", "GBH07")])
```


### Citation:

Stoll & Newbery (to be submitted), Including tree spatial extension in the evaluation of neighbourhood competition effects in Bornean rain forest. Ecol Monogr.

### To do:

- [x] Start preparing appendix (GitHub Repository)
- [x] Make a file list
- [x] Explain function of R-Scripts
- [x] re-work R-Scripts, comments
