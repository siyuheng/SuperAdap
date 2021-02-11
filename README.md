# SuperAdap: The Two-Stage Programming Method in Matched Observational Studies

## Author
Siyu Heng, Hyunseung Kang, Dylan S. Small, and Colin B. Fogarty

## Maintainer
Siyu Heng (Email: <siyuheng@sas.upenn.edu>)

## Description
**SuperAdap** is R package for performing a sensitivity analysis in matched observational studies with the two-stage programming method proposed in the paper "Increasing Power for Observational Studies of Aberrant Response: An Adaptive Approach" by Heng, Kang, Small, and Fogarty. 

Before installing this R package, please ensure that you have installed the following three R packages: **gurobi**, **Matrix**, **mvtnorm**. To install this package in R from GitHub, please run the following commands:

```
install.packages("devtools") 
library(devtools) 
install_github("siyuheng/SuperAdap")
```
## Reference
Heng, S., Kang, H., Small, D. S., and Fogarty, C. B. (2020). Increasing power for observational studies of aberrant response: An adaptive approach. arXiv:1907.06770.
