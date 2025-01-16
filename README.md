Package: PLNMFG
Type: Package
Title: PLNMFG: Pseudo-Label Guided Non-negative Matrix Factorization Model with Graph Constraint for Single-cell Multi-omics Data Clustering
Version: 1.0
Author: Hui Yuan <huiy5533@gmail.com>
Maintainer: Yushan Qiu<yushan.qiu@szu.edu.cn>
Description: This pakege implemented the PLNMFG algorithm that fusion different omics data, which aims to achieve multi-omics cell clustering.
        
### Environment: 
matlab;
### Introduction: 
The algorithm is implemented by matlab.  You need to install the matlab version used in this experiment is 9.11.0.1837725 (R2021b) Update 2.  The experimental data are anno omics data.  The steps are as follows：

###### The possible range of parameter values.  
> eta_range = [0,1,2,3,4,5];
> K_range = [50,100,150];
> delta_range = [0.0001,0.001,0.01,0.1,1];
> beta_range = [0.0001,0.001,0.01,0.1,1];
> r_range = [0,1,2,3,4,5];
> alpha1_range = [0.01,0.1,1,10,100,200,300];
###### Import PLNMFG algorithms from data sources.  
> load('anno.mat'); 
###### Load the model.
> [Y,Q,C,G,U] = PLNMFG(X, gt, option);
##### Load evaluation indicators.
> accuracy = sum(permutedLabels == gt) / length(gt);
> nmi = compute_NMI(gt, permutedLabels);
> ami = AMI(gt, permutedLabels);  
> ari = ARI(gt, max(gt), permutedLabels, max(permutedLabels));


 - If you have any problem, please contact Hui Yuan(m13265355293@163.com) and Yushan Qiu (yushan.qiu@szu.edu.cn)!
