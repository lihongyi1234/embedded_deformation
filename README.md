# C++ implementation of: Embedded Deformation for Shape Manipulation

## What is in this repository?
This repository contains the C++ implementation of embedded deformation (ED) [[1]](#link-to-the-original-paper). The formulation of the different cost functions used in ED are defined with the associated jacobians. In contrasts to [[1]](#link-to-the-original-paper), the optimization uses [Levenberg-Marquardt](http://ceres-solver.org/nnls_solving.html#levenberg-marquardt) from CERES instead of Gauss-Newton.

Several other cost functions are also defined, such as the minimization of model to model similarly as the one from Elastic Fusion (in the journal version).



## Illustration
![example](https://github.com/rFalque/embedded_deformation/raw/master/images/screenshot.png "example of embedded deformation")

## todo
* make the code single thread
* clean up the mesh visualization

## Link to the original paper:
[1] [Embedded Deformation for Shape Manipulation](https://graphics.ethz.ch/~sumnerb/research/embdef/Sumner2007EDF.pdf)
