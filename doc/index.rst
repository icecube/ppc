
.. _ppc-main:

Photon Propagation Code (ppc)
=============================

ppc is a software program that can be used as either a stand-alone command-line executable, or within the IceCube software chain as a module to generate and propagate photons from either in-situ light sources (flasher LEDs, standard candle lasers) or photons emitted or caused by high-energy particles or processes (Cherenkov photons and delayed luminescence).

Files
-----

- gpu/

  ppc for CPU/CUDA

- ocl/

  ppc for OpenCL

- par/

  ppc for std::par

- llh/

  llh/DirectFit

- bfr/

  birefringence simulation code; simulates diffusion patterns in polycrystalline ice

- ice/

  contains selected ice models, configured for use by ppc or llh/DirectFit.


Configuration
-------------

ppc is conifigured via environmental variables, configuration files in the ice directory, and module parameters (or command-line parameters for stand-alone executable).

Environmental variables
+++++++++++++++++++++++

these are set as usual within a shell (with an "export" as necessary). Within a python script these are set with os.putenv(). If the ppc module is loaded with a call to load() these must be set prior to loading the module. If using an import statement instead the environmental variables can be set at any point prior to the Execute() statement.

- system-wide

  these are settings used by the GPU driver that are not directly accessed by ppc

  - **DISPLAY**

    *use: DISPLAY=:0*

  - **COMPUTE**

    *use: COMPUTE=:0*

  - **GPU_DEVICE_ORDINAL**

    *example: GPU_DEVICE_ORDINAL=0,2,3*

    sets which GPUs to use (and makes those not set invisible to the program). In the example above, assuming a system with 4 GPUS 0, 1, 2, and 3, sets to use 3 GPUs 0, 2, and 3, and makes GPU #1 invisible to the program. This re-numbers GPUs: in the example above the GPU numbers that the program sees will be 0, 1, and 2.

  - **CUDA_VISIBLE_DEVICES**

    similar to GPU_DEVICE_ORDINAL, used on systems with CUDA driver

  - **CUDA_CACHE_DISABLE**

    *use: CUDA_CACHE_DISABLE=1*

    force kernel re-compilation each time ppc is loaded/run (for OpenCL version of ppc when run on systems with CUDA driver)

  - **GSL_RNG_SEED**

    *example: GSL_RNG_SEED=$RANDOM*

    sets the gsl random number generator seed, used by llh/DirectFit

- ppc

  - **PPCTABLESDIR**

    *example: PPCTABLESDIR=ice/*

    sets directory where ice and most other configuration files are located. By default uses current directory, ""

  - **NEXTGENDIR**

    *example: NEXTGENDIR=ice/*

    sets directory where new sensor configuration files (om.conf, om.dirs, om.map, om.wv_*.?) are located. By default uses the value of PPCTABLESDIR

  - **ICEMODELDIR**

    *example: ICEMODELDIR=ice/*

    sets directory where ice configuration files (icemodel.par, icemodel.dat, and icemodel.bbl) are located. By default uses the value of PPCTABLESDIR

  - **TILTMODELDIR**

    *example: TILTMODELDIR=ice/*

    sets directory where ice tilt configuration files (tilt.par and tilt.dat or tilt.set and tilt.map) are located. By default uses the value of ICEMODELDIR

  - **PPCHOLEICE**

    *example: PPCHOLEICE=ice/as.dat*

    sets the angular sensitivity file. By default is set to "PPCTABLESDIR/as.dat"

  - **FLWL**

    *FLasher WaveLength*

    *example: FLWL=400*

    overrides the extension of the wv.dat file (containing the wavelength sampling curve). So, in this example the wavelength sampling curve is read in from file wv.400

  - **WFLA**

    *Wavelength of FLAsher*

    *example: WFLA=337*

    overrides the wavelength sampling curve conatined in wv.dat with a single wavelength in nm

  - **FLDR**

    *FLasher DiRection*

    *example: FLDR=-1*

    set FLDR=x+(n-1)*360, where 0<=x<360 and n>0 to simulate n LEDs in a symmetrical n-fold pattern, with first LED centered in the direction x. Negative or unset FLDR simulate a symmetric in azimuth pattern of light. Ignored in particle simulation.

  - **OFLA**

    *"Omit FLAsher"*

    *use: OFLA=0*

    setting this disables the default mode where photons that come back to the flashing DOM are omitted. 

  - **FWID**

    *Flasher beam WIDth*

    *example: FWID=9.7*

    sets the width (in degrees) of the 2d gaussian (von-Mieses-Fisher distribution) that determines the light emission profile of flasher LEDs. Set to -1 to simulate an isotropic emission profile; or set to -2 to simulate a linear (in cos(angle to axis)) cone of light. If greater than 999.0, the accurate lab-measured profile is simulated. Additionally if greater than 1050.0, the value is taken as azimuthal direction to cable (in degrees, value of 1080 is azimuthal direction along x axis), which is placed as a perfectly absorbing cylinder of radius of 2.3 cm, touching the DOM outer surface. The cable orientation specified here is used to block photons at the emission point (before starting the propagation through ice; to specify cable near receiving DOMs use configuration file dx.dat (which normally should match the value specified here for the same DOM). When greater than 999.0 also a perfectly absorbing horizontal harness belt of width of 6.8 cm is simulated, and only 10% of photons exiting the DOM sphere below the equator are retained (unless also using option FTLT).

  - **FZCR**

    *Flasher Z-CorRection*

    *example: FZCR=2.6*

    sets the correction to the LED elevation angle (nominally equal to -0.2 degrees for horizontal and 48.1 degrees for tilted flasher LEDs). As of mid-2025 the best value used was +0.2 degrees for horizonatal and -6.1 degrees for tilted LEDs. Default: 0.

  - **ELEV**

    *flasher led ELEVation angle*

    *example: ELEV=5*

    This sets the flasher LED elevation angle directly (instead of as a correction to the nominal values). This still uses the internal sensor LED coordinates for DOM LEDs, possibly to be changed in the future. So, consider this option experimental for now. Default: not used.

  - **FTLT**

    *Flasher om TiLT*

    *example: FTLT="0 0.6 0.8"*

    This sets the flasher DOM tilt direction vector, which nominally is assumed to lie along the z-axis. If this is used but set to anything other than 3 valid float-point numbers, a DOM tilt value from cx.dat file is used instead. Default: (0, 0, 1).

  - **OVSZ**

    *OVerSiZe factor*

    *example: OVSZ=5*

    overrides the DOM overize scaling factor specified in the cfg.txt file.

  - **BFRA/BFRB**

    *BireFRingence A/B factors*

    *example: BFRA=1.1*

    sets overall (i.e., deflection+diffusion)/deflection-only multiplicative scaling factors for the birefringence effect. Nominal values are 1. These combine (in a product) with the values (if specified) in icemodel.dat configuration file.

  - **HIFL**

    *Hole Ice at/near FLasher DOM*

    *example: HIFL=-20*

    use to only simulate hole ice column(s) near the emitter for "young" photons, when negative (so, younger than 20 ns in the example above), or for older photons (older than then specified value) when positive. This is used in hole ice position fits, to decouple the effects of the hole ice near the emitter vs. near the receivers. Default value of 0 disables this special treatment.

  - **GECO**

    *GEometry COrrection*

    *example: GECO="63 20 0.045 0.17 0.5 0.03"*

    updates x,y position of a specified DOM, if number of given elements is at least 4. The corrections x and y are given in units of DOM radii (16.51 cm). Additional elements, when present, override the values given in the cfg.txt file for the hole ice radius (also, given as a multiple of DOM radii), and effective scattering length in m. In the above example the coordinated of DOM 63,20 are adjusted: x by 0.74 cm and y by 2.8 cm, the hole ice column radius is set to 8.3 cm, and effective scattering length is set to 3 cm. The last two numbers set hole ice properties everywhere in the detector, not just near the specified DOM. This option is mainly used in hole ice fits (to both the DOM positions within the hole and fits to the hole ice column radius and scattering length of ice within the column).

  - **VTHK**

    *Vertical ice layer THicKness*

    *example: VTHK=1*

    Estimate local change in ice layer thickness from the default (of usually 10 m) due to ice layer tilt, when set to 1. Default is 0 (assume the same thickness of ice layers everywhere in the detector regardless of x and y coordinates).

  - **BFRM**

    *BireFRingence scattering correction Method*

    *example: BFRM=1*

    Specifies the method for offsetting the added scattering due to birefringence from the main (Mie) scattering table: 0 (default) is single-value subtraction that keeps the existing wavelength parameterization for overall scattering intact (but is perhaps unphysical); 1 subtracts the birefringence effect from the scattering table values (given at 400 nm) and applies existing wavelength dependence to resulting reduced coefficients - this will requre re-parameterization of the wavelength dependence; and 2: does not subtract the birefringence effect at all, and assumes that the table values of Mie scattering coefficients are are exactly that (to be perhaps used in the future ice fits).

  - **ABSM**

    *ABSorption correction Method*

    *example: ABSM=1*

    Specifies the method for offsetting the intrinsic ice absorption from the dust (Mie) absorption table: 0 (default) assumes the icemodel.dat file specifies the dust (Mie) absorption table; 1 assumes the icemodel.dat file specifies the full absorption at 400 nm, which contains both dust (Mie) and intrinsic ice absorption. Both methods have the same wavelength parameterization and should be equivalent when used with the correctly compiled icemodel.dat file.

  - **NPHO/NPHO_X**

    *example: NPHO_2=512*

    sets the average number of photons to process in a single thread on a GPU. If underscore syntax is used, the number that follows the underscore sets the GPU for which to apply this setting. Setting this to 0 takes that GPU out of use. Default is 1024.

  - **NTHR**

    *Number of THReads*

    *example: NTHR=$[84*1024]*

    sets the number of threads to be run in parralel. This is equivalent to nthr*nblk within the CUDA or OpenCL version of the code and is only available for std::par. It is needed since there is no other means of determining available resouces by standard means available from std::par. If this number is too large, the kernel execution might exceed time limits imposed by the OS (a few seconds) and the program execution will fail. If it is too short the program progress will be slower than optimal. In the example above it is set for an A40 GPU, which has 84 MPs.

  - **HQUO**

    *Hits QUOtient*

    *example: HQUO=4*

    reduces memory footprint of the array of photon hits by the given integer number. This may (rarely) cause a buffer overflow and lost hits (with a printed error), if changed from the default value of 1. This appears to be necessary on gpu001 (either because of the A40 GPU or AMD CPU when used with CUDA). Only applies to the CUDA version of ppc.

  - **XMLT/XMLT_X**

    *example: XMLT_1=4*

    oversubscribes to resouces advertised by the driver by the factor specified. Default is 1 for NVidia and 8 for AMD cards. Only applies to the OpenCL version of ppc.

  - **OCPU/OGPU/OACC**

    *use: OGPU=1*

    use: only CPUs, only GPUs, or only Accelarator cards when specified. Can be combined to specify multiple devices. If not set the program will use every device available to it. In a system with GPUs it is recommended to use OGPU=1. Otherwise, if the driver advertises both GPUs and CPU, the load will be spread equally between both. GPUs, usually being faster, will complete their load quickly and then wait for the CPU device to complete, thus leading to idling GPUs and slower overall execution. Only applies to the OpenCL version of ppc.

  - **BADMP/BADMP_X**

    *example: BADMP=3*

    Specifies the hardware number of a stream multiprocessor (MP) that should not be used. If you are getting warnings "Bad ... MP", that could be due to failing hardware. Try excluding one of the MPs mentioned in the warnings, especially if one is repeating multiple times. This only applies to the CUDA version of ppc, and has only been shown to work on older architectures (such as GTX 295). This can still be used on newer architectures to reduce the load on a GPU (e.g., on a GPU with 20 MPs the load will be reduced from 100% to 19/20=95%), but exclusion of a specific hardware MP might not be guaranteed.

- llh/DirectFit

  - **CYLR**

    *CYLindRical cable*

    *use: CYLR=0/1*

    for cable simulation: (1) simulate straight cylindrical cable, which is faster or (0) curved gaussian-like shape of cable that curves around the DOM and asymptotically approaches the DOM axis above/below the DOM, which is slower (and is the default)

  - **ANGR**

    *ANGular-Restricted*

    *example: ANGR="0 0 1"*

    sets nx ny nz components of the cascade/particle direction. At the same time the angular width of the proposal distribution is set to 0, so the direction is held fixed during iterations. This is overriden if the input file "ini" exists and is successfully read at initialization

  - **FSEP**

    *Flasher SEParation*

    *example: FSEP=1*

    llh sum only includes DOMs that are more than FSEP DOMs away from flasher. Default is 1 (so if DOM 4 is flashing, DOMs 3 and 5 are not used)

  - **SREP**

    *Simulated event REPetitions*

    *example: SREP=10*

    simulate event this many times at each step of the calculation. Default is 1

  - **DREP**

    *Data event REPetitions*

    *example: DREP=250*

    the data file contained averages for this many events. Default is 1. Numbers above 1 are usually used only if there were multiple in-situ light source events taken with the same configuration (e.g. 250 flasher events)

  - **LOOP**

    *number of LOOPs*

    *example: LOOP=1000*

    number of llh steps in a sub-chain. Different search methods might use this number differently. E.g., localized random search has this many evaluations. However, it is repeated 10 times in method 11 (usually used for cascade reconstruction). Default is 1000

  - **NORM**

    *NORMalize*

    *example: NORM=1*

    normalize all waveforms to 1, i.e., only use timing and isolate and exclude the total per-DOM charge information from the likelihood. Default is 0 (disabled).

  - **FAIL**

    *FAIL on warnings*

    *use: FAIL=0/1*

    set to 1 to cause the program to fail on some warnings. Default is 0

  - **FAST**

    *FASTer calculation*

    *use: FAST=0/1*

    1: only use time-integrated charges during simultaneous t0 (start time) and energy unfolding steps. This was shown to produce more stable result, although occasionally somewhat worse llh values. 0: use time-binned charges in parts of the calculation when optimizing t0 and unfolding energy/flasher brightness. Default is 0

  - **UNFM**

    *UNFolding Method*

    *use: UNFM=0/1*

    1: Use slower (but more precise) minimizer (using golden ratio algorithm) for t0 search. Default is 0

  - **MONO**

    *MONOtonize*

    *use: MONO=0/1*

    Set MONO=1 to bin data such that the bins rise and fall monotonically around a single maximum bin. This is done via a modification to the bayesian block algorithm. This works well for flasher data, but you may need to disable this (default behavior) if multiple peaks in data waveforms are important (and the fit needs to reproduce this - such as in tau double bang studies). Default is 0

  - **PEAK**

    *slow/precise PEAK search*

    *use: PEAK=0/1*

    Set PEAK=1 to use a more accurate but slower bin optimization that searches for the best peak bin location in bayesian-block single-peak bin solutions. Otherwise (nominally) the peak bin is assumed to contain the max value bin of the input binning. The captured binned range extends (up to) 500 ns from the front of the peak bin into earlier times, and (up to) 1000 ns from the end of the peak bin into later times. For waveforms with only one bin (e.g., when captured charge is low) this means the bin is extended to cover the full time range of the event (or 5000 ns, whichever is larger), as compared to usually smaller window of only 1500 ns. This can be used with or without MONO=1 option. Default is 0

  - **MLPD**

    *MiLliPeDe*

    *use: MLPD=0/1*

    short for millipede. Enables/disables pattern unfolding: loss profile along the track (0 to reconstruct as a cascade, 1 to reconstruct as a track), or azimuthal flasher light emission profile (1 to enable unfolding into 2 up/down components and 72 azimuthal components spaced out 5 degrees apart; or 0 to use emission profile determined by the FLDR setting). Flasher unfolding has additional options. Below is a full list (+10 options, i.e., 10, 11, 12, 13, and 14 are usually x10 faster for the unfolding stages as most LEDs are seeded with brightness expressed in "number of photon bunches" of e=0.5 from the "ini" file). MLPD=-1 is a "legacy" option, thas is even faster and slightly less precise than MLPD=0 (cascades) or MLPD=10 (flashers).

    ::

       -1:  1 component
       0:   1 component
       1:   2+72 down/up+72 5-deg bins
       2:   2+1 (down/up+led)
       3:   2+1 (back/forward+led)
       4:   1+1 (isotropic+led)
       +10: take e from file "ini", otherwise ini=5

  - **FLSH**

    *FLaSHer*

    *example: FLSH=63,20*

    invokes the flasher mode. Sets the flasher position to the value of the parameter

  - **FDUR**

    *Flasher pulse DURation*

    *example: FDUR=70*

    width of flasher emission pulse in ns assuming rectangular profile. Default is 70 ns

  - **QSAT**

    *Q (charge) SATuration value*

    *example: QSAT=500*

    maximum integrated charge per DOM to accept that DOM into the calculation. Default is 500

  - **CNUM**

    *Cos-bin NUMber*

    *example: CNUM=40*

    number of cos(arrival angle wrt. PMT axis) bins. Default is 1

  - **LSIG**

    *Llh SIGma*

    *example: LSIG=0.05*

    value of the sigma/model error to be used in the likelihood evaluation. 0 reverts to likelihood containing only Poisson terms. Default it 0.1 (i.e., 10%)

  - **FLOR**

    *FLasher board ORientation*

    *example: FLOR=1*

    1: tilt the flasherboard in a direction consistent with the DOM tilt (only when MLPD=1). Default is 0

  - **APFR**

    *AfterPulse FRaction*

    *example: APFR=0.01*

    add afterpulses to simulated photons. Expected value is around 1 - 2%. Only a single template is implemented, only accurate up to ~4 us.

- inv (the code used to fit RDEs and to unfold angular sensitivity curve)

  - **SREP**
  - **DREP**

    these have the same meaning as when used with llh and described in the previous section

  - **IGEO**

    *Input GEOmetry*

    *example: IGEO=ice/geo-f2k*

    sets the geometry file

  - **IEFF**

    *Initial EFFiciencies*

    *example: IEFF=ice/eff-f2k*

    eff-f2k file used in the simulation

  - **IORI**

    *Input ORIginal efficiencies*

    *example: IORI=ice/eff-f2k.ori*

    sets the file specifying nominal RDE values (for use with XMAX and XSIG parameters described below)

  - **IANG**

    *Input ANGular sensitivity file*

    *example: IANG=ice/as.dat*

    sets the angular sensitivity file used in the simulation

  - **XINI**

    *X INItial value*

    *example: XINI=xini*

    sets the file containing an appoximation to the unfolding result (usually a result from the previous interation)

  - **XMAX**

    *X MAXimum deviation*

    *example: XMAX=1.5*

    sets the hard limits on RDEs around values contained in the IORI file. In the example above the limits are [1/1.5; 1.5] for a DOM with a nominal RDE value of 1. 0 disables the hard limits. Default is 0

  - **ESCL**

    *Efficiency SCaLing*

    *example: ESCL=1.1*

    scales DOM efficiencies by amount specified here before running the fit. This might be needed to counteract the bias due to the llh model.

  - **XSIG**

    *X-SIGma*

    *example: XSIG=0.01*

    adds regularization around the RDE values specified in the IORI file with width specified. Default is 0.1

Configuration files
+++++++++++++++++++

- ice (set by PPCTABLESDIR)

  - **as.dat**

    Definition of the angular sensitivity of a DOM. There are two possible variations in the format:

      - 1st number is greater than 0:

        the rest of the file contains coefficients of the polinomial expansion of the angular sensitivity curve vs. cos(photon arrival angle wrt. PMT axis). The first number is the maximum value reached by this curve (and is applied as a cut regardless of whether it's the actual high point of the polinomial curve on the interval -1,1 or not). Numbers lower than 1 accelerate the calculation (since fewer photons need to be simulated)

      - 1st number is 0:

        This defines the "surface sensitivity" option. Second number in the file defines sensitive area: cos(angle to photon hit point on surface from center wrt. PMT axis) must be greater than this number to accept the photon

  - **cfg.txt**

    main configuration file. See example below for explanation.

    ::

      # ppc configuration file: follow strict order below
      5     # over-R: DOM radius "oversize" scaling factor
      1.0   # overall DOM efficiency correction
      0.35  # f_SL (scattering function) 0=HG; 1=SAM
      0.9   # g=<cos(theta)> (scattering function)

      130.0 # azimuth of major anisotropy axis (deg)
      0.0   # magnitude of major anisotropy coefficient k1
      0.0   # magnitude of minor anisotropy coefficient k2

      0.5   # hole ice radius in units of [DOM radius]
      0.03  # hole ice effective scattering length [m]
      100   # hole ice absorption length [m]
      0.35  # hole ice 0=HG; 1=SAM
      0.9   # hole ice g=<cos(theta)>

      0.6   # magnitude of major anisotropy coefficient k1
      -0.3  # magnitude of minor anisotropy coefficient k2
      -0.3  # magnitude of minor anisotropy coefficient kz

      0.0   # scaling for old absorption anisotropy

      0.076795 # p1, sigma along flow
      544284.5 # p2
      2.229494 # p3
      0.002624 # p4
      0.077381 # p1, sigma perpendicular to flow
      1547618. # p2
      2.449589 # p3
      0.002505 # p4
      0.000995 # p1, mean deflection towards the flow
      0.248264 # p2
      2.354436 # p3
      1.680717 # p4

    The block of 12 birefringence parameters p1/p2/p3/p4=A/B/C/D=alpha/beta/gamma/delta is used as described in the New Anisotropy Paper (formulae 1,2) and in Cryosphere (formulae 23,24). See list of references at the end of this page for links to references. Blocks other than the first one can be optionally omitted (disabling anisotropy and hole ice parts of the calculation)

  - **cx.dat**

    DOM tilt map, each line contains: String#, OM#, nx, ny, nz, uncertainty (degrees). nx, ny, nz are components of the tilt vector that is defined as opposite of PMT axis direction

  - **dx.dat**

    Cable position map, each line contains: String#, OM#, azimuth direction to cable (degrees), uncertainty (degrees).

  - **eff-f2k**

    RDE (relative DOM efficiency) map, each line contains: String#, OM#, RDE, Type. If no entry RDE=1, Type=0 are assumed. DOMs that use corrected  wavelength acceptance from file wv.rde (for high-QE DOMs) have Type=1. It is possible to specify high-QE DOMs with Type=0 and simply a higher RDE value (nominally 1.35), of with an RDE value near 1 and Type=1. The acceptance correction curve parametrized in wv.rde file nears a value of 1.35 for wavelengths near 400 nm. RDE values taken from the GCD frame are matched with Type=0. If a corrected wavelength dependence is desired, GCD values need to be overridden by having this file (and wv.rde) present in the ice configuration directory

  - **geo-f2k**

    Geometry map, each line contains: DOM ID, Mainboard ID, x, y, z, String#, OM#. This file is necessary for running ppc from command-line. When present and running as an icetray module, will override the values from GCD

  - **str-f2k**

    String geometry map; specifies x,y coordinates of centers of drilled holes (presumed to coincide with the center of the bubbly column of hole ice running the length of the string). Each line contains: String#, x, y. This file is necessary when geo-f2k file specifies precise locations of DOMs, to include position within the drilled holes, i.e., relative to the bubbly (central) column of the hole ice. When not present, the average of x,y coordinates of DOMs 1-60 on each string is used instead.

  - **hvs-f2k**

    High-voltage map, each line contains: String#, OM#, high voltage. Used only to specify that the DOM is on when HV>0. This file overrides the map of "ON" DOMs from GCD when present in the ice directory.

  - **icemodel.bbl**

    parametrization of air bubble contribution to scattering. Has 3 values: b, d1, d2. The parametrized contribution is b*(d1-d)*(d2-d) for d that specifies a shallower depth than both d1 and d2. The contribution is 0 otherwise (i.e. for deeper locations)

  - **icemodel.dat**

    main ice properties table: depth of the center of the layer, be(400), adust(400), delta tau (as defined in section 4 of the SPICE paper). All layers must be of equal width, and there must be at least 1 layer defined in the file. If there is only one layer, the layer width is set to 100 m, otherwise it is calculated from the layer depth spacing. If the file icemodel.par contains 6 parameters, then the absorption coefficient is calculated as adust(400)=(D*[3rd element in a line]+E)*400^-kappa.

    This file may contain 2 additional optional columns, containing the anisotropy coefficients k1 and k2. Ice layers defined with lines containing k1 and k2 will use these anisotropy coefficients instead of those specified in file cfg.txt

    Finally, two more columns can be present, containing the birefringence "strength" parameters BFRA/BFRB. These set overall (i.e., deflection+diffusion) and deflection-only multiplicative scaling factors for the birefringence effect.

  - **icemodel.par**

    file with 4 parameters of the icemodel: alpha, kappa, A, B (as defined in section 4 of the SPICE paper). Each parameter is followed by its measurement uncertainty, which is ignored by the program. The older models (older than SPICE Lea or WHAM) have 6 parameters: alpha, kappa, A, B, D, E.

  - **rnd.txt**

    table of random number multipliers for the multiply-with-carry random number generator used by the parallelized kernel. Can have one or more elements per line, but only the first number is used (this is to make is copmatible with older formats of this file)

  - **tilt.dat**

    Describes ice layer tilt, each line contains: center depth of layer, and several depth corrections for locations specified in file tilt.par

  - **tilt.par**

    Containes a map of tabulated ice layer tilt locations, each line has: string number, and a relative distance along the gradient tilt direction (225 degrees SW)

  - **tilt.map**

    Describes new (2d) ice layer tilt, each line contains: center depth of layer, followed by depth corrections for locations configured in file tilt.set. When files tilt.map and tilt.set are present and configure a valid tilt model at initialization, files tilt.dat and tilt.par are ignored.

  - **tilt.set**

    Containes a grid configuration for tabulated ice layer tilt locations, each line (of two) has: azimuthal direction of the grid axis, coordinate of the first element, step size, and number of elements along the axis. To define a hexagonal region, the last element configures the exclusion region (where the SE and NW corners are cut to create the two extra sides of the hexagon). Only two directions for the two axes are currently supported: 9.3 and 129.3 degrees (these are chosen to align with IceCube detector geometry). When files tilt.map and tilt.set are present and configure a valid tilt model at initialization, files tilt.dat and tilt.par are ignored.

  - **wv.dat**

    parametrization of wavelength-tabulated DOM acceptance (calculated from qe_dom2007a table of efficiency.h file of photonics), convolved with input spectrum. Each line contains: normalized integrated acceptance, and wavelength in nm.

  - **wv.rde**

    parametrization of the correction to the wavelength acceptance curve to be used for high-QE DOMs. Each line has: wavelength in nm, and correction factor (ratio of high-QE to nominal)

  - **om.conf**

    main OM (Optical Module) configuration file. See example below for explanation.

    ::

       # name  module  area    beta    Rr      Rz      num     dir     cable
       Default -1      1       0.49    0.1651  0.1651  1       180 0   90
       DEgg    120     1       0.5     0.150   0.267   2       180 0
                                                                 0 0

    Area is an overall efficiency scaling parameter for a given sensor type. Parameter beta specifies the angular sensitivity function of a single PMT within any given sensor; values between -1 and 1 specify PMTs with approximately spherical cathodes (including flat with beta=1); values less than -1 (e.g., -2) specify cylindrical sensitive area. Parameters Rr and Rz specify extensions of the spheroid that approximates the sensor geometry in the horizontal and vertical directions (spheroid is assumed symmetrical around the z-axis). If Rz is negative, the sensor is assumed cylindrical instead (with the radius Rr and length of -2*Rz). Num is the number of PMTs within the sensor, and dir specify zenith and azimuth of each PMT. Cable, when specified, is the azimuthal direction to cable in the sensor's coordinate system. The sensor orientation itself is then determined from this number combined with the azimuthal direction to cable as given in file dx.dat.

  - **om.dirs**

    contains a uniformly distributed set of directions on a sphere; needed to calculate normalization factors (maximum sensitivity/projected area and direction-averaged area vs. maximum projected area) of sensors at initialization.

  - **om.map**

    OM (Optical Module) type map, each line contains: String#, OM#, OM type. At the moment all PMTs within the same sensor are assumed to be of same (PMT) type, which is specified as before in file eff-f2k.

  - **om.wv_[om_type].[pmt_type]**

    parametrization of wavelength-tabulated direction-averaged sensor effective area. Each line contains: direction-averaged sensor effective area in cm^2, and wavelength in nm. OM type in filename extension is taken from om.conf and PMT type (all PMTs within the same sensor are at the moment assumed to be of the same type) is taken from file eff-f2k.

- llh/DirectFit additional configuration/input files, to be placed in the "current" directory

  - **as**

    this has the same format as as.dat in the ice directory. llh needs to be able to apply the angular sensitivity within its code when fitting for the DOM tilt or cable position. When using file "as" make sure to apply a uniform/flat sensitivity curve in file as.dat (e.g., by having it contain 2 numbers: 0.68 and 0.68)

  - **zs**

    contains the grid of search directions, used when fitting the DOM tilt (which is performed if this file is found). Each line contains: a unique identifying number, nx, ny, nz. This file can be generated with program "ico" in llh subdirectory

  - **cx**

    this has the same format as cx.dat in the ice directory. Make sure that only one of "cx", "cx.dat" is available at run time. If fitting for DOM tilt iteratively with llh/DirectFit, make sure that only "cx" is available.

  - **cs**

    contains the set of azimuthal positions of cable to test used when fitting for the cable position (which is performed if this file is found). Each line contains: a unique identifying number, and azimuth angle in degrees.

  - **dx**

    this has the same format as dx.dat in the ice directory. Make sure that only one of "dx", "dx.dat" is available at run time. If fitting for cable position iteratively with llh/DirectFit, make sure that only "dx" is available.

  - **bad**

    contains String#, OM# of DOMs that are to be considered bad in the fit. If this file is found the DOMs in it are excluded from the fit, and the no-hit contribution to the log likelihood sum is taken into account.

  - **ert**

    contains String#, OM#, ti, tf that define the "DOM errata" list containing time intervals of bad data, which are not to be used. May define more than one interval [ti; tf) for each DOM.

  - **set**

    contains String#, OM# of the "good DOM list". When present and read, only uses DOMs that are on this list.

  - **dat**

    main data file. Each line contains: String#, OM#, time in ns, and average charge in p.e.s. The data is internally rebinned in 25 ns bins before applying the bayesian deblocking method to merge bins. If an event spans over more than 5000 ns then to avoid resizing the fixed 200 bin internal buffers the bin size is increased. It is recommended to trim events to keep then at 5000 ns or less in length but throwing away late pulses and coinsident events before or after the main event. Coincident events should be cleaned away anyway with, e.g., topological trigger. Longer events such as muons crossing the entire detector should of course not be shortened just to fit into 5000 ns, but only to remove afterpulses and coincident events.

  - **sim**

    Uses the same format as the "dat" file above. Overrides the binned charge values while re-using the binning calculated/optimized for the "dat" file. This is useful for re-simulations of the data, where it might be desired to keep the same binning as used for the actual data.

  - **ini**

    Contains cascade/flasher parameters (to be used as best fit, or as initial approximation, or to facilitate iterations passing the solution between separate runs of llh). It may contain one or more lines, ordered as listed below:

    1) x, y, z (meters), zenith, azimuth (degrees), energy (GeV)/flasher brightness (bunches), time (ns), scattering and sbsorption scaling coefficients (last two unsupported, set both to 1.0). This line needs 5 or more elements to be accepted (some values at the end may be omittes, like the scaling coefficients)
    2) sequence numbers representing the unfolded pattern. The number of elements must match that expected by llh (usually defined by the geometry of the event) exactly, and the elements should sum up to 1. This line may be left empty if it is not needed but the following lines are.
    3) estimates of the proposal distribution parameters: rr (correlation between position and direction), dr (spacial width in m), da (angular width in degrees), di (intended for use with scattering and absorption scaling coeffifients, so should be left as 0), and optionally, lx (threshold value of likelihood, used in the ABC method for calculating uncertainties)

Command-line parameters
+++++++++++++++++++++++

- ppc

  - no parameters

    Prints a summary of available tables, and an error if something is missing. If all necessary tables are found, also prints a summary of the available GPU devices within your system. These are numbered starting with 0 and must be specified with a [gpu] parameter in the examples below.

  - one parameter "-"
  - optionally "-" [x] [y]

    Print out the table of ice parameters (IceCube coordinate z of the center of the ice layer, absorption coefficient, effective scattering *Mie* coefficient, and 4pi-integrated birefringence *BFR* scattering) for wavelength w in [nm] (if set with WFLA=[w]) at the IceCube coordinates x and y in [m] (or 0, 0 if not specified). The parameters are computed using formulae of section 4 of the SPICE paper.

  - one parameter "_"

    Print out a table containing ice layer tilt surfaces

  - one parameter "="

    expects x, y, z of a point in IceCube coordinates. Outputs the input string appended by z-shift (due to tilt) to and local ice layer width (nominally 10 m but can vary if using VTHK=1) at the x,y position of the nominal ice given by the icemodel.dat table. This is used to compile the table matching DOMs to ice layers:

    ::

       #!/bin/sh

       ppc=ppc
       ice=ices/spice_ftp-v3

       awk 'BEGIN {zoff=1948.07} $6>0 && $7<=60 {print $3, $4, $5+zoff, $6, $7}' geo-f2k | PPCTABLESDIR=$ice $ppc = |
       awk 'BEGIN {zoff=1948.07; low=zoff-2798.47} {z=$3-$6; n=1+int((z-low)/10+0.5); print $4, $5, n}' > strx-ftp

  - one integer parameter [gpu]

    Process particle simulation in f2k format from stdin. The muons must have been processed by mmc with the "-recc" option, which prints out all muon segments individually as "amu" particles. An additional (to the f2k spec) TR particle "ph" can be specified, which inserts a photon (or integer number of photons, exactly specified and inserted as is, without sampling a Poisson distribution around given number or scaling for oversize). Here is an example of f2k input to ppc:

    ::

      #!/bin/awk -f

      BEGIN {
        print "V 2000.1.2"
        print "TBEGIN ? ? ?"

        srand(1);
        for(i=0; i<100; i++){
          x=(2*rand()-1)*500 # meters
          y=(2*rand()-1)*500 # meters
          z=(2*rand()-1)*500 # meters
          zenith=rand()*180  # degrees
          azimuth=rand()*360 # degrees
          l=500              # length, m
          energy=1.e5        # GeV
          t=0                # ns
          N=1    # number of photons

          print "EM 1 1 1970 0 0 0"
          print "TR 1 0 e    ", x, y, z, zenith, azimuth, 0, energy, t
          print "TR 1 0 amu  ", x, y, z, zenith, azimuth, l, energy, t
          print "TR 1 0 hadr ", x, y, z, zenith, azimuth, 0, energy, t
          print "TR 1 0 ph   ", x, y, z, zenith, azimuth, 0, N, t
          print "EE"
        }
        print "TEND ? ? ?"
        print "END"
      }


  - 4 parameters: [str] [dom] [num] [gpu]

    Simulate [num] photons emitted by a flasher or a standard candle at the position [str],[om]. Please note the following rules:

    - positive [str] simulates horizontal flashers, negative [-str] simulates tilted flashers,
    - str=0 and om=1,2 simulates standard candles 1 and 2,
    - you must set WFLA=337 before simulating the standard candles,
    - if the wv.dat file contains the flasher wavelength profile, WFLA=405 should be omitted,
    - if [num] is specified as x*y, x photons are simulated y times (with y trailing empty lines).

- llh/DirectFit

  - 1 parameter [method]

    - 0, 1 (same as 0): calculate llh. -1 additionally prints out the simulated hit data for the best solution
    - 10: applies localized random seach after calculating initial guess
    - 16: ABC (Approximate Bayesian Calculation) to estimate uncertainties



Description of output
---------------------

- ppc

  command-line ppc reports hits with "HIT" lines:

  HIT String# OM# time(ns) wavelength(nm) p_theta p_phi d_theta d_phi

  p_theta and p_phi specify direction of the photon at the point where it impacts the DOM

  d_theta and d_phi specify direction from the DOM center to the point of photon impact



- llh

  \*,?,|,-,+,0,1,... llh x y z zenith azimuth n t s a

  first element is a special character or an integer indicating the progress of the program as it goes through iterations. The second element, llh, is the saturated log likelihood (a measure of the goodness of fit), that indicates how well the simulation matches data (lower values are better). Coordinates x,y,z are in meters, zenith and azimuth in degrees. The next element, n, is either reconstructed particle deposited energy in GeV or light source brightness in photon bunches. This number maybe a constant factor off (effects like SPE mean, cable shadow, etc.) One way to figure out this factor is to reconstruct a few simulated events of known energy (i.e., calibrate the output of llh/DirectFit with the proper IceCube simulation). Next is t, the t0 (in ns) of the event. Finally, s and a are the scaling ice corrections. These are currently not used and are left at 1 each.

  If negative method is used as run-time parameter, the best match between data and simulation will be printed out in the following format:

  String# OM# bin_size charge_data(p.e.) charge_simulation(p.e.)

  internally llh/DirectFit bins the data nominally in 25 ns bins. It might be necessary to increase the bin size to a larger number if the overal length of the event is larger than 5 us (200 of the 25 ns bins). The bin size is printed out on stderr at initialization as "lbin". Then llh applies the Bayesian blocks procedure to merge some of the data bins. The number of initial bins contained in such a merged bin is indicated in the line above as bin_size. For each String#,OM# the bins are printed out in their time order, so it should be possible to infer the time structure (waveform) of the detected charge, albeit only with precision limited by the variable-size bins. This can be used to plot data and best simultion for visual inspection of the quality of fit.


Random notes on code structure
------------------------------

the point of next scatter is calculated by solving the following equation:

Exp( - integral_0^(distance) (scattering coefficient at x) dx ) = (random number)

since we have ice organized in layers of constant optical properties the integral reduces to a sum, and calculating the distance to next scatter is as simple as solving a linear equation with a couple boundary checks.


the point of absorption is calculated by solving the following equation:

Exp( - integral_0^(distance) (absorption coefficient at x) dx ) = (random number)

since we have ice organized in layers of constant optical properties the integral reduces to a sum, and calculating the distance to next scatter is as simple as solving a linear equation with a couple boundary checks. One complication compared to scattering is that the sum is done over multiple segments because of intermediate scatterings. So the code keeps subtracting the integral evaluated between successive scatters from the -log(random number) until it drops below zero. When that happens the particle does not make it to the next scatter point and the point of absorption is calculated instead.

on any segment the direction is fixed and absorption coefficient is modified according to the anisotropy model. It should be easy to do the same to the scattering coefficient, if necessary.

it is possible to make ppc create a re-producible sequence of photons by uncommenting "#define DTMN" within the ppc core (files gpu/ini.cxx, ocl/ppc.cxx, par/ini.cxx} or adding -DDTMN to the compiler options. This will make the CUDA and OpenCL versions of the code slower (by ~ 15%) and will only result in the same sequence when the runs are repeated on the same GPU.

Specific code details requested in tickets
------------------------------------------

Cable shadow
++++++++++++

Cable shadow is implemented in ppc as a folloowing approximation: the photon landing coordinates (on a DOM) and final direction are used to "backtrack" the photon to sek whether it could have intersected with cable positioned next to the DOM before landing on the surface of the DOM. The cable shape is implemented as a vertical cylinder with a radius of 23 mm that is touching the surface of the DOM (for the nominally oriented DOM, at the equator). An implementation of the cable as a curve touching the DOM surface but asymptotically aligning with the DOM axis above and below the DOM, also exists within the llh/DirectFit code. This code can be switched on by setting CYLR=0. The location of the cable wrt. the DOM is specified via angle to the cable in the configuration file dx.dat.

The cable shadow code is implemented outside of the propagate kernel (which runs on the GPU), thus executing on the CPU side, during the post-processing of the hits (the code is in file f2k.cxx).

Hole Ice
++++++++

The hole ice is implemented to describe the following physical construction precisely: A vertical infinite cylinder column of constant ice properties, going through the center of the string (for each string). The ice properties are specified in the cfg.txt file in the optional block as described in the configuration files section above. The configuration describes the radius of the cylindrical column, scattering and absorption coefficients, and shape parameters of the scattering function: f_SL and g. If all of the DOMs of the string have exactly the same x and y coordinates, the hole ice column is simulated as concentric with the DOMs. In order to simulate the situation where the hole ice is to a side of the DOM, the DOM coordinates need to be adjusted. Keep in mind that this will in turn modify the average x and y of the string (i.e. of the center of the string), so the coordinates of the rest of the DOMs need to be adjusted in the opposite direction by a small amount (1/60th for a nominal IceCube string). Since this implementation is currently just a placeholder, still awaiting detailed calibration of the hole ice properties, a more verstile configuration has not yet been implemented. It may turn out that the future configurations fully implementing hole ice can be fully specified with the existing scheme, or that modification may be required.

DOM tilt
++++++++

DOM tilt is implemented by assuming a tilt in the DOM axis (i.e. deviation from the vertical) during application of the angular DOM sensitivity. This is done for either the "effective" angular DOM sensitivity that only depends on the direction of the photon (and then the angle on which the angular sensitivity depends is calculated wrt. the DOM tilt axis rather than the vertical), or the new DARD-style DOM sensitivity which accepts photons a certain distance down from the DOM equator (i.e., effectively simulating a sensitive survace of the PMT). The tilt directions of DOMs are specifies in the configuration file cx.dat.

The DOM tilt code is implemented outside of the propagate kernel (which runs on the GPU), thus executing on the CPU side, during the post-processing of the hits (the code is in file f2k.cxx).

DOM Oversize
++++++++++++

In order to accelerate calculation, the DOMs can be (and normally are) oversized. The latest ice models were fitted with the oversize factor of 5, and much of the nominal IceCube simulation also uses an oversize factor of 5. Simplifying the actually ipmlemented geometry (explained in the following paragraph) a little, this means the DOM radius is increased by the factor of 5, increasing the DOM cross-sectional area presented to the oncoming photons by a square of the oversize factor, i.e., by a factor of 25. This, in turn, means that a factor 25 times fewer photons need to be generated and propagated, thus accelerating the simulation by a factor of 25.

A number of variations of the oversizing geometries and simulation strategies were studied some time after the feature was introduced, and the best option was settled on during the development of the SPICE Mie ice model (and used mostly unchanged since then). The main issue with scaling the DOM as a perfect sphere was that the DOM was occupying much more space, 125 time more (a cube of the oversize factor). This space was then unavailable to the scattering and absorption processes, significantly changing the timing distributions (as well as time-integrated charge). Thus, another approach, the "pancake" scheme, was implemented instead. In the pancake construction the DOM dimensions are scaled only in the directions perpendicular to the traveling photon, while maintaining the nominal dimension in the direction along the photon travel. This maintains the factor 25 scaling in the DOM cross section presented to the photon, while also reducing the volume occupied by the oversized DOM from a factor 125 to only 25 times the nominal volume. As the photon scatters and changes direction, so does the pancake rotate so that the area presented to the photon is always 25 times the nominal. Such changes in simulated DOM geometries, as well as the larger dimensions of the DOMs (compared to nominal) do lead to some deviations in the timing distributions, and for oversize factor of 16 (used in the development of ice models up to and including SPICE Lea) lead to about 1-2 ns distortion of the timing distribution at 17 m, and 3-9 ns distortion at 125 m (as measured in the positions of the leading edge and peak in the clear deep ice). SPICE 3.x ice models used the entire volume of IceCube detector in the ice model fits, including the more dense part DeepCore. It was thought that the oversize factor 16 was too big to adequately approximate the physics of the denser parts of the array in cleaner ice, so a smaller factor of 5 (also matching the factor used in most of the simluation produced at the time) was settled on and used in SPICE 3.x models.

Additional feature in the strategy of the oversize implementation is to continue propagating the photon after it hits the DOM. This is disabled when the oversize factor is set to 1, and the photon is stopped and disappears once it hits the surface of the DOM. Such a strategy is thought to be more correct, as in the equivalent nominal-size DOM treatement, while propgating a bunch of 25 photons that would otherwise hit the oversized DOM, only one of them will hit the nominal-size DOM, while 24 will continue unimpeded. This strategy was also included into the timing distribution distortion numbers stated in the previous paragraph.

RDE implementation
++++++++++++++++++

The file wv.rde (ratio of high-QE to nominal vs. wavelength) will be phased out.
New approach taken for the configuration of IceCube extensions is to specify the effective area vs. wavelength for each type of sensor PMT individually.


References:
-----------

SPICE models
++++++++++++

`Optical properties of deep glacial ice at the South Pole (AMANDA ice paper) <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2005JD006687>`_

`Measurement of South Pole ice transparency with the IceCube LED calibration system (SPICE Paper): arXiv:1301.5361 <https://arxiv.org/abs/1301.5361>`_

`Evidence of optical anisotropy of the South Pole ice (Ice Anisotropy Paper): arXiv:1309.7010 (pp. 17-20) <https://arxiv.org/abs/1309.7010>`_

`Light diffusion in birefringent polycrystals and the IceCube ice anisotropy (Birefringence Paper): arXiv:1908.07608 <https://arxiv.org/abs/1908.07608>`_

`A novel microstructure-based model to explain the IceCube ice anisotropy (New Anisotropy Paper): arXiv:2107.08692 <https://arxiv.org/abs/2107.08692>`_

`A calibration study of local ice and optical sensor properties in IceCube (Hole Ice Paper): arXiv:2107.10435 <https://arxiv.org/abs/2107.10435>`_

`In-situ estimation of ice crystal properties at the South Pole using LED calibration data from the IceCube Neutrino Observatory (Cryosphere) <https://tc.copernicus.org/preprints/tc-2022-174/>`_

`A model independent parametrization of the optical properties of the refrozen IceCube drill holes <https://arxiv.org/abs/2307.15298>`_

`An improved mapping of ice layer undulations for the IceCube Neutrino Observatory <https://arxiv.org/abs/2307.13951>`_

llh/DirectFit
+++++++++++++

`Event reconstruction in IceCube based on direct event re-simulation (DirectFit paper): arXiv:1309.7010 (pp. 21-24) <https://arxiv.org/abs/1309.7010>`_

`Likelihood description for comparing data with simulation of limited statistics (Likelihood Paper): arXiv:1304.0735 <https://arxiv.org/abs/1304.0735>`_

`Likelihood description for comparing data to simulation of limited statistics (LLH ICRC Paper) <https://inspirehep.net/literature/1413017>`_

PPC
+++

`Photon tracking with GPUs in IceCube, Nuclear Inst. and Methods in Physics Research, A, Volume 725, pp. 141-143. <https://www.sciencedirect.com/science/article/pii/S0168900212015173>`_

`Photon Propagation with GPUs in IceCube, Proceedins to GPUHEP2014, DESY-PROC-2014-05, pp. 217-220 <https://inspirehep.net/literature/1386642>`_

Miscellaneous
+++++++++++++

`Older ppc information pages on my homepage <https://icecube.wisc.edu/~dima/work/WISC/ppc/>`_ and `readme file <https://icecube.wisc.edu/~dima/work/WISC/ppc/readme.html>`_

`AMANDA file format definition (used as input for command-line ppc in particle mode) <https://www-zeuthen.desy.de/~steffenp/f2000/>`_

`Muon Monte Carlo (MMC), a java program that can be used to process muons for use with command-line ppc <https://icecube.wisc.edu/~dima/work/MUONPR/>`_
