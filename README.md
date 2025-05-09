# MOR-adiabatic-mode
A MATLAB code for solving HREs in the modeling of 3D underwater acoustic propagation

This repository hosts a MATLAB toolkit implementing an adiabatic mode framework for 3D underwater acoustic propagation in longitudinally invariant waveguides. The numerical solution employs modal decomposition with coefficients efficiently computed through a Model Order Reduction (MOR) technique. For complete theoretical details and numerical validation, we refer readers to [JSV 2024](https://doi.org/10.1016/j.jsv.2024.118617) by He, et al. This package is a supplementary material to reproduce the results presented in the [JSV paper](https://doi.org/10.1016/j.jsv.2024.118617).

# INSTALLATION and EXECUTION

-- Download the package

-- Run MATLAB codes

# ASA wedge and cosine hill cases

Exam1_HRESolver_TransSymtr_ASAwedge.m

Exam2_HRESolver_TransSymtr_CosHill.m

Note：

Folders 'ASA_wedge_kj/' and 'CosHill_kj/' contain: 1) the horizontal wavenumbers k_rj, 2) eigenfunctions at the receiver depth, and 3)  mode excitation at the source depth for the wedge and canyon examples, respectively. Benchmark reference solution is archived in the associated folders. Execution of the MATLAB scripts Exam1.m and Exam2.m will regenerate all numerical results underlying Figures 7-13 in the [JSV paper](https://doi.org/10.1016/j.jsv.2024.118617).

# Gaussian canyon case

Exam3_HRESolver_TransSymtr_Canyon.m

Note：

This demonstration reproduces Figs. 8 and 9 presented in [JASA 2021](https://pubs.aip.org/asa/jasa/article/150/2/1140/615453/A-three-dimensional-finite-difference-model-for) by Liu, et al, which is the 3D propagation over an underwater canyon. This case has also been demonstrated previously in [JASA 2019](https://pubs.aip.org/asa/jasa/article/146/3/2050/995175/Split-step-Pade-solver-for-three-dimensional) and [Appl. Sci 2020](https://www.mdpi.com/2076-3417/10/7/2393). The same parameters and geometry were adopted in this code as in these papers. 

Folder 'Canyon_kj/' contains 1) the horizontal wavenumbers k_rj, 2) eigenfunctions at the receiver depth, and 3)  mode excitation at the source depth.

# Internal wave case

Exam4_HRESolver_TransSymtr_internalwaves.m

Note：

This numerical demonstration simulates three-dimensional underwater acoustic propagation through internal wave fields in a stratified marine environment. The experimental configuration replicates [SWARM 95](https://pubs.aip.org/asa/jasa/article/117/2/613/541579/Measurement-and-modeling-of-three-dimensional). A 70 m depth waveguide is considered in this code, where the sound speed profile is consist of three layers: surface isovelocity  layer with a sound speed of 1530 m/s, the thermocline layer with a gradient of 3.33 1/s, and bottm isovelocity layer with a sound speed of 1480 m/s. The thermocline layer initially starts at 15 m depth, with a thickness of 15 m. An internal wave train moves at 0.65 m/s at the upper and lower interfaces of the thermocline layer, with each wave having an amplitude of 10 m and a wavelength of 70 m. The train contains 10 circles of waves. The seabed parameters can be found in 'param.m'. In line 58 of 'Exam4.m', the varaint 'dy' is the whole time sequence when the internal waves occur. Adjusting the index of 'dy(1)' in line 62 can produce the sound field at different time while the internal wave train is moving.

Function 'ModeSolver.m' provides the solver for the eigenequation associated with waterborne modes. This solver can handle the situation where a multi-layered medium is considered. 

Script 'param.m' provides the input to 'ModeSolver.m', which is the protocal where the environment parameters are given. Comments are provided following each parameter, which is easy to read and understand for others.


