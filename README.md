# MOR-adiabatic-mode
A MATLAB code for solving HREs in the modeling of 3D underwater acoustic propagation

This repository hosts a MATLAB toolkit implementing an adiabatic mode framework for 3D underwater acoustic propagation in longitudinally invariant waveguides. The numerical solution employs modal decomposition with coefficients efficiently computed through a Model Order Reduction (MOR) technique. For complete theoretical details and numerical validation, we refer readers to [JSV 2024](https://doi.org/10.1016/j.jsv.2024.118617) by He, et al. This package is a supplementary material to reproduce the results presented in the [JSV paper](https://doi.org/10.1016/j.jsv.2024.118617).

# INSTALLATION and EXECUTION

-- Download the package

-- Run MATLAB codes

Exam1_HRESolver_TransSymtr_ASAwedge.m

Exam2_HRESolver_TransSymtr_CosHill.m

# Note：

Folders 'ASA_wedge_kj/' and 'CosHill_kj/' contain the horizontal wavenumbers k_rj, eigenfunctions at the receiver depth, and mode excitation at the source depth for the wedge and canyon examples, respectively, with the reference solution stored in the associated folders. Running Exam1 and Exam2 direclty reproduces the results presented in Figs. 7-13 of the [JSV paper](https://doi.org/10.1016/j.jsv.2024.118617).

Exam3_HRESolver_TransSymtr_Canyon.m

# Note：

This demonstration reproduces Figs. 8 and 9 presented in [JASA 2021](https://pubs.aip.org/asa/jasa/article/150/2/1140/615453/A-three-dimensional-finite-difference-model-for) by Liu, et al, which is the 3D propagation over an underwater canyon. This case has also been demonstrated previously in [JASA 2019](https://pubs.aip.org/asa/jasa/article/146/3/2050/995175/Split-step-Pade-solver-for-three-dimensional) and [Appl. Sci 2020](https://www.mdpi.com/2076-3417/10/7/2393). The same parameters and geometry were adopted in this code as in these papers. 

Folder 'Canyon_kj/' contains the horizontal wavenumbers k_rj, eigenfunctions at the receiver depth, and mode excitation at the source depth.

Exam4_HRESolver_TransSymtr_internalwaves.m

# Note：


