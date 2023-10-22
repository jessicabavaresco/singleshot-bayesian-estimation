## Code to accompany: *[title](https://arxiv.org/abs/xxxx.xxxxx)*
#### Jessica Bavaresco, Patryk Lipka-Bartosik, Pavel Sekatski, and Mohammad Mehboudi

This is a repository for the code used in the article "*title*, Jessica Bavaresco, Patryk Lipka-Bartosik, Pavel Sekatski, and Mohammad Mehboudi, [arXiv:xxxx.xxxxx [quant-ph]](https://arxiv.org/abs/xxxx.xxxxx)".

All code is written in MATLAB and requires:
- [Yalmip](https://yalmip.github.io) - a free MATLAB toolbox for rapid prototyping of optimization problems
- [MOSEK](https://www.mosek.com) - a software package for solving mathematical optimization problems (under the free personal academic license)
- [QETLAB](http://www.qetlab.com/) - a free MATLAB toolbox for quantum entanglement theory

This repository consists of the following:

#### CODE

- [SDP_scoreoptimization.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/SDP_scoreoptimization.m):
**SDP** that solves the optimization in **Eq. (17)**. That is, that given a set of operators { $X(\hat{\theta}_i)$ }, maximizes or minimizes (depending on the reward/cost function) the score $\mathcal{S}$ over testers $T=$ { $T_i$ }.

- [script_phaseestimation.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/script_phaseestimation.m):
Script for the metrology problem presented in ***Sec. V.A. Example 1: Paradigmatic example – Local phase estimation*** It evaluates methods 1, 2, and 3 for the parameters of this example, outlined in the paper.

- [script_thermometry.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/script_thermometry.m):
Script for the metrology problem presented in ***Sec. V.B. Example 2: Non-unitary evolution – Thermometry***. It evaluates methods 1, 2, and 3 for the parameters of this example, outlined in the paper.

- [script_su2estimation.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/script_su2estimation.m):
Script for the metrology problem presented in ***Sec. V.C. Example 3: Multi-parameter estimation – SU(2) gates***. It evaluates methods 1, 2, and 3 for the parameters of this example, outlined in the paper.

#### DATA FILES

- [Data/data_phaseestimation.mat](https://github.com/jessicabavaresco//singleshot-bayesian-estimation/blob/main/Data/data_phaseestimation.mat):
MATLAB data file generated with [script_phaseestimation.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/script_phaseestimation.m) and plotted in **Fig. 2** of the paper. It contains the  variables
  - score_M1, score_M2, score_M3: vector containing the best score found applying methods 1, 2, and 3, respectively. Each line corresponds to one value of $N_o\in$ { $2,..,10$ }.
  - T_M1, T_M2, T_M3: cells containing the best testers found applying methods 1, 2, and 3, respectively. Each cell entry corresponds to a tensor that contains the best tester for a value of $N_o\in$ { $2,..,10$ }. The first two indeces of each tensor correspond to a matrix defining a tester element and the third index marks the label $i\in$ { $1,..,N_0$ } of the tester element. 

- [Data/data_thermometry.mat](https://github.com/jessicabavaresco//singleshot-bayesian-estimation/blob/main/Data/data_thermometry.mat):
MATLAB data file generated with [script_thermometry.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/script_thermometry.m) and plotted in **Figs. 3-5** of the paper. It contains the variables
   - score_M1, score_M2, score_M3: matrix containing the best score found applying methods 1, 2, and 3, respectively. Each line corresponds to a value of $N_o\in$ { $2,..,10$ } and each column corresponds to a time step.
  - T_M1, T_M2, T_M3: cells containing the best testers found applying methods 1, 2, and 3, respectively. Cells: each line corresponds to a value of $N_o\in$ { $2,..,10$ } and each column corresponds to a time step; each entry contains the a tensor with the best tester found in that scenario. The first two indeces of each tensor correspond to a matrix defining a tester element and the third index marks the label $i\in$ { $1,..,N_0$ } of the tester element.
  - score_M1_PPT, score_M2_PPT, score_M3_PPT, T_M1_PPT, T_M2_PPT, T_M3_PPT: analogous variables for the no-entanglement approximation of the problem.

- [Data/data_su2estimation.mat](https://github.com/jessicabavaresco//singleshot-bayesian-estimation/blob/main/Data/data_su2estimation.mat):
MATLAB data file generated with [script_su2estimation.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/script_su2estimation.m) and plotted in **Fig. 6** of the paper. It contains the variables
  - score_M1, score_M2, score_M3: vector containing the best score found applying methods 1, 2, and 3, respectively. Each line corresponds to one value of $N_o\in$ { $2,..,10$ }.
  - T_M1, T_M2, T_M3: cells containing the best testers found applying methods 1, 2, and 3, respectively. Each cell entry corresponds to a tensor that contains the best tester for a value of $N_o\in$ { $2,..,10$ }. The first two indeces of each tensor correspond to a matrix defining a tester element and the third index marks the label $i\in$ { $1,..,N_0$ } of the tester element. 

