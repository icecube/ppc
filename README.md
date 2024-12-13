# ppc
photon propagation code (ppc) and selected IceCube ice models

- ppc sources are contained in two directories that are fully equivalent to each other in functionality:  
   gpu: for CPU (one thread) or CUDA  
   ocl: for OpenCL  
   par: for std::par version of the code

- llh: llh/DirectFit. Compiles against either version of ppc (gpu or ocl).  
   use this to fit ice models to flasher data or reconstruct cascade/muon IceCube events

- bfr: BireFRingence code  
   calculates birefringence patterns for various crystal ice configurations

- ice: contains bare-minimum configurations for selected IceCube ice models (to be used with ppc)  
   spice_bfr-v2: South Pole ICE based on BireFRingence with absorption anisotropy model

- doc: [documentation](doc/index.rst)


<a href="https://zenodo.org/doi/10.5281/zenodo.10410725"><img src="https://zenodo.org/badge/527689072.svg" alt="DOI"></a>
