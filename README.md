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
You need to install Yalmip, SeDuMi, and COMPLeib to run the main program.
- Yalmip can be downloaded from [https://yalmip.github.io/](https://yalmip.github.io/)
- SeDuMi can be downloaded from [https://sedumi.ie.lehigh.edu/](https://sedumi.ie.lehigh.edu/)
- COMPLeib can be downloaded from [http://www.complib.de/](http://www.complib.de/)
