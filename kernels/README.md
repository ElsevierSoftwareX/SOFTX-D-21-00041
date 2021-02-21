# Elastodynamic Kernels

to create kernels for spectral boundary integral method solving the elastodynamic wave equation

![](.images/convolution_kernels.png)

*Note that H22 is independent of nu*

## Usage

For generating H00 H01 H11 kernels

	./laplace_inverse.py <nu> [pstress]
		
For generating H22 kernel

	./laplace_inverse.py h22



## Requirements

* `python-mpmath`
* `python-numpy`
* `python-matplotlib`
* `python-tk`

## Notes

Please be aware that the generation of the kernels is (in python) a very time consuming process. However, this is only required if you are considering a material for which the kernel file does not yet exist. Hence, this is rarely used.
