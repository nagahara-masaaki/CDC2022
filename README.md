# Matlab programs
We provide Matlab programs to check the numerical example of
- M. Nagahara and N. Sebe, Reduced-order controller design via iterative greedy LMI, IEEE CDC 2022

## Main program
The main programs are `reduced_order_controller_CDC2022_stability.m` (for Section IV.A) and `reduced_order_controller_CDC2022_hinf_bisection.m` (for Section IV.B).
## Subprograms
You need to place the following m-files in the same folder of the main programs.
- `reduced_order_controller_CDC2022_stability.m` (for Section IV.A)
- `reduced_order_controller_CDC2022_hinf_bisection.m` (for Section IV.B)
- `Proj_LMI.m`
- `Proj_rank.m`
- `retrieve_eliminated.m`
- `s_sparse_operator.m`
## Requirements
You need to install Yalmip, SeDuMi, and COMPLeib to run the main program.
- Yalmip can be downloaded from [https://yalmip.github.io/](https://yalmip.github.io/)
- SeDuMi can be downloaded from [https://sedumi.ie.lehigh.edu/](https://sedumi.ie.lehigh.edu/)
- COMPLeib can be downloaded from [http://www.complib.de/](http://www.complib.de/)
