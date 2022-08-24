# Matlab programs
We provide Matlab programs to check the numerical examples of model reduction (MR) shown in
- M. Nagahara, Y. Iwai, and N. Sebe, Projection onto the Set of Rank-constrained Structured Matrices for Reduced-order Controller Design, submitted to Algorithms, 2022.

## Main program
The main programs are `reduced_order_controller_stability.m` (for Section 4.1) and `reduced_order_controller_hinf_bisection.m` (for Section 4.2).
## Subprograms
You need to place the following m-files in the same folder of the main programs.
- `reduced_order_controller_stability.m` (for Section 4.1)
- `reduced_order_controller_hinf_bisection.m` (for Section 4.2)
- `Proj_LMI.m`
- `Proj_rank.m`
- `retrieve_eliminated.m`
- `s_sparse_operator.m`
## Requirements
To run `hinfstruct` for nonsmooth H-infinity synthesis, you need [Robust Control Toolbox](https://jp.mathworks.com/products/robust.html)
You also need to install Yalmip, SDPT3, SeDuMi, and COMPLeib to run the main program.
- Yalmip can be downloaded from [https://yalmip.github.io/](https://yalmip.github.io/)
- SDPT3 can be downloaded from [https://github.com/SQLP/SDPT3](https://github.com/SQLP/SDPT3)
- SeDuMi can be downloaded from [https://sedumi.ie.lehigh.edu/](https://sedumi.ie.lehigh.edu/)
- COMPLeib can be downloaded from [http://www.complib.de/](http://www.complib.de/)
## Results of stability tests for 95 bemchmark models
The results of stability tests for 95 benchmark models are shown in the following exel file:
- `stability_data.xlsx`
