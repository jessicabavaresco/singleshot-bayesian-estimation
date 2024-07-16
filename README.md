## Code to accompany: *[Designing optimal protocols in Bayesian quantum parameter estimation with higher-order operations](https://arxiv.org/abs/2311.01513)*
#### Jessica Bavaresco, Patryk Lipka-Bartosik, Pavel Sekatski, and Mohammad Mehboudi

This is a repository for the code used in the article "*Designing optimal protocols in Bayesian quantum parameter estimation with higher-order operations*, Jessica Bavaresco, Patryk Lipka-Bartosik, Pavel Sekatski, and Mohammad Mehboudi, [*Phys. Rev. Research* **6**, 023305 (2024)](https://doi.org/10.1103/PhysRevResearch.6.023305), [arXiv:2311.01513 [quant-ph]](https://arxiv.org/abs/2311.01513)".

All code is written in MATLAB and requires:
- [Yalmip](https://yalmip.github.io) - a free MATLAB toolbox for rapid prototyping of optimization problems
- [SeDuMi](https://github.com/sqlp/sedumi) - a linear/quadratic/semidefinite solver for Matlab and Octave
- [QETLAB](http://www.qetlab.com/) - a free MATLAB toolbox for quantum entanglement theory

This repository consists of the following:

#### CODE

- [SDP_scoreoptimization.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/SDP_scoreoptimization.m):
**SDP** that solves the optimization in **Eq. (18)**. That is, that given a set of operators { $X(\hat{\theta}_i)$ }, maximizes or minimizes (depending on the reward/cost function) the score $\mathcal{S}$ over testers $T=$ { $T_i$ }.

- [tester2realization.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/tester2realization.m):
Function that computes a quantum state and measurement that compose a quantum realization of a given tester, according to Eqs. (12) and (13) of the paper.

- [script_phaseestimation.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/script_phaseestimation.m):
Script for the metrology problem presented in ***Sec. V.A. Example 1: Paradigmatic example – Local phase estimation*** It evaluates methods 1, 2, and 3 for the parameters of this example, outlined in the paper.

- [script_thermometry.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/script_thermometry.m):
Script for the metrology problem presented in ***Sec. V.B. Example 2: Non-unitary evolution – Thermometry***. It evaluates methods 1, 2, and 3 for the parameters of this example, outlined in the paper.

- [script_su2estimation.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/script_su2estimation.m):
Script for the metrology problem presented in ***Sec. V.C. Example 3: Multi-parameter estimation – SU(2) gates***. It evaluates methods 1, 2, and 3 for the parameters of this example, outlined in the paper.

#### DATA

- [data_phaseestimation.mat](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/data_phaseestimation.mat):
Output data obtained from [script_phaseestimation.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/script_phaseestimation.m) containing the optimal scores, testers, states and measurements found by our methods for the metrology problem presented in ***Sec. V.A. Example 1: Paradigmatic example – Local phase estimation***.

- [data_thermometry_scores&testers.mat](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/data_thermometry_scores&testers.mat):
Output data obtained from [script_thermometry.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/script_thermometry.m) containing the optimal scores and testers found by our methods for the metrology problem presented in ***Sec. V.B. Example 2: Non-unitary evolution – Thermometry***.

- [data_thermometry_states.mat](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/data_thermometry_states.mat):
Output data obtained from [script_thermometry.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/script_thermometry.m) containing the optimal states found by our methods for the metrology problem presented in ***Sec. V.B. Example 2: Non-unitary evolution – Thermometry***.

- [data_thermometry_measurements.mat](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/data_thermometry_measurements.mat):
Output data obtained from [script_thermometry.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/script_thermometry.m) containing the measurements found by our methods for the metrology problem presented in ***Sec. V.B. Example 2: Non-unitary evolution – Thermometry***.

- [data_su2estimation.mat](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/data_su2estimation.mat):
Output data obtained from [script_su2estimation.m](https://github.com/jessicabavaresco/singleshot-bayesian-estimation/blob/main/script_su2estimation.m) containing the optimal scores, testers, states and measurements found by our methods for the metrology problem presented in ***Sec. V.C. Example 3: Multi-parameter estimation – SU(2) gates***.


