# Robust topology optimization for structures under bounded random loads and material uncertainties

Thank you for your interest in this project!

This repo contains part of the codes for the current research by Song Bai and Zhan Kang from Dept. of Engineering Mechanics, Dalian University of Technology, Dalian, China.

<https://doi.org/10.1016/j.compstruc.2021.106569>

Please read the notes below for additional pre-requisites for the implementation of the codes.

NOTE:

1. The MATLAB version of the Method of Moving Asymptotes (MMA) optimizer code ought to be obtained by contacting Prof. Krister Svanberg.

2. The subroutine "ParforProgressbar.m" is the open-source code developed by Frerk Saxen to monitor parallel computation progress. It needs to be obtained independently via the link: <https://www.mathworks.com/matlabcentral/fileexchange/71436-parfor-progress-monitor-progress-bar-v4>.

3. The subroutine "nsumk.m" is used to calculate the multi-index of PCE base functions. It needs to be obtained via <https://www.mathworks.com/matlabcentral/fileexchange/28340-nsumk>.

4. The subroutine "nwspgr.m" is employed to generate Smolyak sparse grid. It needs to be obtained via <http://www.sparse-grids.de/>.

5. The subroutine "setprod.m" is adopted to calculate the tensor-product quadrature nodes' indices. It is to be obtained via <https://www.mathworks.com/matlabcentral/fileexchange/5898-setprod>. However, since the sparse grid scheme is employed instead of the tensor product scheme, "setprod.m" is not actually used.

## Topology evolution animations of the numerical examples

* Numerical example 1: Bridge beam
![Example 1](./images/topology_evo_example_1.gif)

* Numerical example 2: Planar structure
![Example 2](./images/topology_evo_example_2.gif)

* Numerical example 3: Cuboid structure
![Example 3](./images/topology_evo_example_3.gif)
