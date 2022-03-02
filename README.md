# Matlab programs
We provide Matlab programs to check the numerical example of
- M. Nagahara and N. Sebe, Reduced-order controller design via iterative greedy LMI, IEEE CDC 2022

## Main program
The main program is `reduced_order_controller_CDC2022.m`
## Subprograms
You need to place the following m-files in the same folder of the main program.
- `controller_realization.m`
- `Proj_LMI.m`
- `Proj_rank.m`
- `retrieve_eliminated.m`
- `s_sparse_operator.m`
## Requirements
You need to install Yalmip and SeDuMi to run the main program.
- Yalmip can be downloaded from [https://yalmip.github.io/](https://yalmip.github.io/)
- SeDuMi can be downloaded from [https://sedumi.ie.lehigh.edu/](https://sedumi.ie.lehigh.edu/)
