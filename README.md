# ppc
photon propagation code (ppc) and selected IceCube ice models

- ppc sources are contained in two directories that are fully equivalent to each other in functionality:  
   gpu: for CPU (one thread) or CUDA  
   ocl: for OpenCL

- llh: llh/DirectFit. Compiles against either version of ppc (gpu or ocl).  
   use this to fit ice models to flasher data or reconstruct cascade/muon IceCube events

- bfr: BireFRingence code  
   calculates birefringence patterns for various crystal ice configurations

- ice: contains bare-minimum configurations for selected IceCube ice models (to be used with ppc)  
   spice_bfr-v2: South Pole ICE based on BireFRingence with absorption anisotropy model

- doc: [documentation](doc/index.rst)
