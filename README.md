# Group Structure Estimation for Non-linear Panel Data 



This repository provides the implementation of the method described in (Yu et al., 2021). 

## Usage
The file `MembershipEstimation.R` contains a function `membership.est` to estimate the membership using the proposed approach. The file `NumGrpsEstimation.R` contains a function `G.est` to estimate the number of groups using the proposed heuristic. The aforementioned paper contains detailed descriptions of the methods. 

## Simulations
We provide an example code `Model2_MembershipEstimation.R` for examining the performance of the proposed group structure estimation method in terms of perfect match (the proportion of times that the exact group assignment is found) and average match (average percentage of correct classification, see the aforementioned paper for more details). The data is generated from Model 2 in (Yu et al., 2021). 

The example code `Model2_NumGrpsEstimation.R` is provided for examining the performance of the proposed heuristic ( eq (2) in Yu et al., 2021) for the number of groups estimation. The data is generated from Model 2 in (Yu et al., 2021).


## Repository Authors
- Lu Yu — Department of Statistical Sciences, University of Toronto 

- Jiaying Gu — Department of Economics, University of Toronto

- Stanislav Volgushev — Department of Statistical Sciences, University of Toronto 


## Reference
[Yu, L., Gu, J, and Volgushev, S. (2021) Group structure estimation for non-linear panel data -- a general approach.]
