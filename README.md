# KEPLER
Kernelized Saliency-Based Person Re-Identification Through Multiple Metric Learning - Transactions on Image Processing (2015)

This package provides you an updated version of the MATLAB code for following paper: 
Martinel, N., Micheloni, C., & Foresti, G. L. (2015). Kernelized Saliency-Based Person Re-Identification Through Multiple Metric Learning. IEEE Transactions on Image Processing, 24(12), 5645–5658. http://doi.org/10.1109/TIP.2015.2487048

## USAGE:

### Simple Demo
```MATLAB
 startup;
```
Initialize all the directiories which are needed to run the main algorithm

```MATLAB
results = main();
```
Reproduces a single experiment using the VIPeR dataset. Results, like the CMC and the normalized Area Under Curve are stored in the results structure.

### Settings
If you want to play with the approach parameters, please refer to the 
```MATLAB
init_parameters.m
```
file.


## ADDITIONAL TOOLBOX:
With this package some additional libraries used in the method are also provided. Note that the algorithm works with the given libraries versions and it's not guaranteed to work with newer or older ones.

    vlfeat: A. Vedaldi and B. Fulkerson - VLFeat: An Open and Portable Library of Computer Vision Algorithms (http://www.vlfeat.org/)
    KISSME: M. Köstinger, M. Hirzer, P. Wohlhart, P. M. Roth, H. Bischof - Large Scale Metric Learning from Equivalence Constraints (https://lrs.icg.tugraz.at/research/kissme/)

## COMPILE:
Please note that some libraries contain mex-files that needs to be compiled for your machine. We provide a limited set of binary within the package. We have not yet provided a script to compile all the dependences.

## CITATION:
If you use the code contained in this package we appreciate if you'll cite our work. 
BIBTEX:
@article{Martinel2015c,
author = {Martinel, Niki and Micheloni, Christian and Foresti, Gian Luca},
doi = {10.1109/TIP.2015.2487048},
issn = {1057-7149},
journal = {IEEE Transactions on Image Processing},
month = {dec},
number = {12},
pages = {5645--5658},
title = {{Kernelized Saliency-Based Person Re-Identification Through Multiple Metric Learning}},
url = {http://ieeexplore.ieee.org/lpdocs/epic03/wrapper.htm?arnumber=7289409},
volume = {24},
year = {2015}
}
