Scripts to generate optimized Runge-Kutta methods for a given spatial discretization. This uses `RK-Opt` (https://github.com/ketch/RK-Opt). Note that this repo contains scripts that I have written to use `RK-Opt`; for any questions about `RK-Opt` itself, consult its documentation (https://rk-opt.readthedocs.io/en/latest/index.html) and/or contact its authors.

Cloning
-------

This repository contains `RK-Opt` as a submodule. Either clone using

```
git clone --recurse-submodules https://github.com/debog/OptimizedRK.git
```
or 
```
git clone https://github.com/debog/OptimizedRK.git
git submodule init
```

Running
-------

Run the script `ConsFD_1stDeriv`. The directory `example_output` contains an example solution of a 2nd order, 5-stage explicit RK method that is optimized for the 5th-order WENO scheme. One can change the parameters in `ConsFD_1stDeriv` to try various orders, number of stages, etc.
