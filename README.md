# WaterMeanMeanForce

Constrained and restrained path integral molecular dynamics (PIMD) for the computation of the quantum potential of mean force (PMF) of water clusters.
Currently tailored to the dimer system.

Tested with Julia 1.5.


## Installation

```
pkg> add https://github.com/0/WaterMeanMeanForce.jl.git
```

In order to run the driver scripts in `bin/`, you will also need to
```
pkg> add ArgParse
```

### Application project

If you're working with a clone of this repository, you can use the basic application project in `bin/`, which already has `WaterMeanMeanForce` and `ArgParse` as dependencies.
From the repository root, run
```
julia --project=bin
```
and then
```
pkg> dev .
```
to create `bin/Manifest.toml` with a development version of `WaterMeanMeanForce`.


## Examples

To run the following examples, you should set the project (e.g. using `--project` or `JULIA_PROJECT`) to a Julia project that has the prerequisites installed.

* `julia bin/harmonic_exact_slope.jl -T 10,50,100`
* `julia bin/run_equilibration.jl --model q-SPC/Fw --integrator OBABO -T 10 -R 0.75 --num-links 16 --time-step 1e-4 --centroid-friction 10 --equil-duration 10 --umbrella-k 1000 > equil.dat`
* `julia bin/run_production.jl --model q-SPC/Fw --integrator c-BAOAB -T 10 -R 0.75 --num-links 16 --time-step 1e-4 --centroid-friction 10 --equil-duration 10 --duration 10`


## References

The binning method is from: Ambegaokar, V., & Troyer, M. (2010). Estimating errors reliably in Monte Carlo simulations of the Ehrenfest model. American Journal of Physics, 78(2), 150-157. [doi:10.1119/1.3247985](https://aapt.scitation.org/doi/abs/10.1119/1.3247985), [arXiv:0906.0943](https://arxiv.org/abs/0906.0943).

The q-SPC/Fw water model is from: Paesani, F., Zhang, W., Case, D. A., Cheatham III, T. E., & Voth, G. A. (2006). An accurate and simple quantum model for liquid water. The Journal of Chemical Physics, 125(18), 184507. [doi:10.1063/1.2386157](https://aip.scitation.org/doi/abs/10.1063/1.2386157).

The OBABO (PILE) integrator is from: Ceriotti, M., Parrinello, M., Markland, T. E., & Manolopoulos, D. E. (2010). Efficient stochastic thermostatting of path integral molecular dynamics. The Journal of Chemical Physics, 133(12), 124104. [doi:10.1063/1.3489925](https://aip.scitation.org/doi/abs/10.1063/1.3489925), [arXiv:1009.1045](https://arxiv.org/abs/1009.1045).

The BAOAB integrator is from: Liu, J., Li, D., & Liu, X. (2016). A simple and accurate algorithm for path integral molecular dynamics with the Langevin thermostat. The Journal of Chemical Physics, 145(2), 024103. [doi:10.1063/1.4954990](https://aip.scitation.org/doi/abs/10.1063/1.4954990), [arXiv:1611.06331](https://arxiv.org/abs/1611.06331).

The post-quantization restraint (path integral umbrella sampling) method is from: Bishop, K. P., & Roy, P.-N. (2018). Quantum mechanical free energy profiles with post-quantization restraints: Binding free energy of the water dimer over a broad range of temperatures. The Journal of Chemical Physics, 148(10), 102303. [doi:10.1063/1.4986915](https://aip.scitation.org/doi/abs/10.1063/1.4986915).

The weighted histogram analysis method (WHAM) is from: Kumar, S., Rosenberg, J. M., Bouzida, D., Swendsen, R. H., & Kollman, P. A. (1992). The weighted histogram analysis method for free‐energy calculations on biomolecules. I. The method. Journal of Computational Chemistry, 13(8), 1011-1021. [doi:10.1002/jcc.540130812](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540130812).


## Acknowledgements

Thanks to Kevin Bishop for helping to verify this implementation!


## License

Provided under the terms of the MIT license.
See `LICENSE` for more information.
