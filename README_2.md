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

Run the script `CompactFD_3rdDeriv`. The directory `Compact_output` contains an example solution of a 3rd order, 4-stage explicit RK method that is optimized for the 8th-order TDCCS and TDCNCS scheme. 
