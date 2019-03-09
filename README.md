# itpfit
This is a code for computing Inflationary Three-Point Functions Involving Tensors (ITPFIT) 
in canonical single field inflationary models.

---
## Using this code

You can use this code provided that you cite the following papers in your publications.

	1. [Numerical evaluation of the three-point scalar-tensor cross-correlations and the tensor bi-spectrum](https://doi.org/10.1088/1475-7516/2013/12/037)
	2. [Examining the consistency relations describing the three-point functions involving tensors](https://doi.org/10.1088/1475-7516/2014/10/021)

---
## About this code

1. This code consists of two files, viz. main.f90 (the main program) and 
   qpws.f90 (a module).
2. The file main.f90 calculates the scalar-tensor three-point functions and 
   the tensor bi-spectrum for an inflationary model involving a single, 
   canonical scalar field for an arbitrary triangular configuration of the 
   wavevectors k1, k2 and k3. In the process, it also evaluates the scalar 
   as well as the tensor power spectra, which are required to arrive at the 
   corresponding non-Gaussianity parameters.
3. All the information about the potential is contained in the module. As an 
   example (as mentioned above), we have provided a module for the quadratic 
   potential with a step (qpws.f90). In order to use the code for another 
   canonical single field model, the module qpws.f90 needs to be modified 
   suitably.

4. The initial value of scale factor (ai) is arrived at by demanding that the 
   pivot scale leaves the Hubble radius at 50 efolds before end of inflation. 
   In case, a model does not lead to the end of inflation within, say, 60-70 
   e-folds, then ai may have to be defined suitably.

---
## Compiling and running the code

In order to run the code, in a terminal, simply type the following:

gfortran/ifort -mcmodel=medium qpws.f90 main.f90
./a.out

---
## Support

To get support, please open a new issue on: 
[https://github.com/v-sreenath/itpfit](https://github.com/v-sreenath/itpfit)
