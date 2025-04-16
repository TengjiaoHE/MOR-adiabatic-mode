# MOR-adiabatic-mode
A MATLAB code for solving HREs in the modeling of 3D underwater acoustic propagation
This repository contains MATLAB package for the adiabatic solution of 3D underwater acoustics propagation in longitudinally invariant environments. The solution is represented in the form of the modal decomposition, and the respective modal amplitudes are obtained using model order reduce (MOR) technique. More details refer to the [JSV paper](https://doi.org/10.1016/j.jsv.2024.118617). This package is a supplementary material to reproduce the results presented in the [JSV paper](https://doi.org/10.1016/j.jsv.2024.118617).

INSTALLATION and EXECUTION

-- Download the package

-- Run MATLAB codes

Exam1_HRESolver_TransSymtr_ASAwedge.m

Exam2_HRESolver_TransSymtr_CosHill.m

Exam3_HRESolver_TransSymtr_Canyon.m

Exam4_HRESolver_TransSymtr_internalwaves.m


Note. Folders 'ASA_wedge_kj/' and 'canyon/' contain wavenumbers k_j and eigenfunction values at z=z_r for the second (wedge) and third (canyon) examples, respectively. Alson the data with reference solutions is stored in these folders.
