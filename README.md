- [Report new Bug Issue](https://github.com/RSE-Cambridge/qe-gpu/issues/new)
- [Guidelines for Contributing](CONTRIBUTING.md)
- [Pull Request Template](.github/PULL_REQUEST_TEMPLATE.md)
- [Project License](License)


## GPU-accelerated Quantum ESPRESSO (QE-GPU)

This is an open-source custom version of Quantum ESPRESSO with embedded GPU
support based on CUDA FORTRAN. This product has been made possible thanks to
the effort of the [NVIDIA](http://www.nvidia.com/page/home.html) HPC Software
and Benchmarks Group. This version is maintained by
[Filippo Spiga](https://github.com/fspiga) and hosted on the University of
Cambridge Research Software Engineering GitHub page. Partial support was
provided by [E4 Computer Engineering SpA](https://www.e4company.com/en/) via
the European PRACE Pre-Commercial Procurement project (Phase 3). To contribute
please refer to the guidelines in [CONTRIBUTING.md](CONTRIBUTING.md)


### Requirements

The [PGI](http://www.pgroup.com/products/community.htm) compiler version 17.4 
or above is required to use QE-GPU. It containes CUDA SDK 8.0 and pre-built 
Open MPI for parallel execution (check the
[PGI Instalation Guide](http://www.pgroup.com/doc/pgiinstall174.pdf) how to 
install it). **no other compilers are supported**

You need data-centre grade NVIDIA TESLA Kepler (K20, K40, K80) or NVIDIA TESLA
Pascal (P100) compute GPUs. No other cards are supported. NVIDIA TESLA P100 is
strongly recommend for its memory capacity and performance.

This version of QE-GPU runs **exclusively** in parallel, Open MPI is required
and also Intel MKL.


### Installation

To compile QE-GPU, copy a `make.inc` template from "install/" directory into the main directory and run make.

These make.inc templates are available:
* `make.inc_x86-64` to compile on any x86-64 machines with NVIDIA GPU (`GPU_ARCH={35, 60}`)
* `make.inc_x86-64_CPU-only` to compile on x86-64 without GPU support
* `make.inc_CRAY_PizDaint` to compile on Piz Daint at CSCS, CRAY XC30 with P100 PCI GPU support (`GPU_ARCH=60`)
* `make.inc_POWER_DAVIDE` to compile on PRACE D.A.V.I.D.E. at E4, POWER8 system with GPU support (`GPU_ARCH=60`)
* `make.inc_POWER_SUMMITDEV` to compile on ORNL early access system in preparation of next OLCF's next big supercomputer, SUMMIT (`GPU_ARCH=60`)

By invoking _make_ alone a list of acceptable targets will be displayed. Binaries go in "bin/". Read comments in the `make.inc` templates to customize it further based on your ebvironment and where math libraries are located. The architectures/environments supported are x86-64, POWER and CRAY.

The QE-GPU package has been reduced in size to the minimum essential. For more
information, please refer to the general documentation provided with the full
Quantum ESPRESSO suite or visit the official web site
[http://www.quantum-espresso.org/](http://www.quantum-espresso.org/)


### Citation

Zenodo: [![DOI](https://zenodo.org/badge/80047177.svg)](https://zenodo.org/badge/latestdoi/80047177)



### License

All the material included in this distribution is free software; you can
redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

These programs are distributed in the hope that they will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 675 Mass
Ave, Cambridge, MA 02139, USA.
