# FBM_ADT
This is the source code of the paper "Reliability modeling and statistical analysis of accelerated degradation process with memory effects and unit-to-unit variability"

If you find the code useful, please give a star :)

## Ref

Chen, Shi-Shun, Xiao-Yang Li, and Wen-Rui Xie. "Reliability modeling and statistical analysis of accelerated degradation process with memory effects and unit-to-unit variability." Applied Mathematical Modelling 138 (2025): 115788.

## Highlights

• A non-Markovian accelerated degradation model with memory effect is developed.

• Unit-to-unit variability is considered in the accelerated degradation model.

• Ignoring unit-to-unit variability leads to biased estimation of the memory effect.

• A statistical analysis method is proposed via expectation maximization algorithm.

• The proposed statistical analysis method gives a more accurate estimation.

Paper link: https://www.sciencedirect.com/science/article/pii/S0307904X24005419

The folder "code_simulation" contains the source code of the simulation study.

In the folder, `main_sim_FBM_utuv.m` is the proposed model M0; `main_sim_FBM.m` is the model M1; `main_sim_BM_utuv.m` is the model M2; `main_sim_BM.m` is the model M3.

## **Attention!** **Corrigendum**

We have identified an error in the figures 8, 10 and 11 in our article, which is related to the illustration of the deterministic trend of the degradation prediction. This error was accidentally introduced due to the picture modification in the revision process, and affected only the pictures. The numerical results in the tables and conclusions remain accurate.

We sincerely apologize for this oversight. The corrected figures are provided below.

![image](https://github.com/dirge1/FBM_ADT/blob/main/corrected_figures/fig.8.png)
Fig. 8. Prediction for the deterministic degradation trend under all stress levels.

![image](https://github.com/dirge1/FBM_ADT/blob/main/corrected_figures/fig.10.png)
Fig. 10. Prediction for the cross-validation 1: (a) deterministic degradation trend

![image](https://github.com/dirge1/FBM_ADT/blob/main/corrected_figures/fig.11.png)
Fig. 11. Prediction for the cross-validation 2: (a) deterministic degradation trend
