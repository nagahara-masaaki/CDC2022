# Matlab programs for reduced-order controller design
We provide Matlab programs to check the numerical examples of model reduction (MR) shown in Sections 4.1 and 4.2 of
- M. Nagahara, Y. Iwai, and N. Sebe, Projection onto the Set of Rank-constrained Structured Matrices for Reduced-order Controller Design, submitted to Algorithms, 2022.

## Matlab programs
You need to place the following m-files in the same folder of the main programs.
- `reduced_order_controller_stability.m` (main program for Sections 4.1 and 4.2)
- `reduced_order_controller_hinf_bisection.m` (main program for Section 4.3)
- `Proj_LMI.m` (projection operator for the LMIs)
- `Proj_rank.m` (projection operator for the rank condition)
- `controller_realization.m` (controller realization from X and Y)
- `retrieve_eliminated.m` (used in `controller_realization.m`)
- `s_sparse_operator.m` (s-sparse operator for truncation)
## Requirements
To run `hinfstruct` for nonsmooth H-infinity synthesis, you need [Robust Control Toolbox](https://jp.mathworks.com/products/robust.html).
You also need to install Yalmip, SDPT3, SeDuMi, and COMPLeib to run the main program.
- Yalmip can be downloaded from [https://yalmip.github.io/](https://yalmip.github.io/)
- SDPT3 can be downloaded from [https://github.com/SQLP/SDPT3](https://github.com/SQLP/SDPT3)
- SeDuMi can be downloaded from [https://sedumi.ie.lehigh.edu/](https://sedumi.ie.lehigh.edu/)
- COMPLeib can be downloaded from [http://www.complib.de/](http://www.complib.de/)
## Results of stability tests for 95 bemchmark models
The results of stability tests for 95 benchmark models are shown in the following exel file:
- `stability_data.xlsx`
