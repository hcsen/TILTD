# TILTD


## Setup

### mice (Linux)

You will need suitable toolchains installed.
On Mahuika `module load MATLAB/{version}` is enough.

```
wget https://naif.jpl.nasa.gov/pub/naif/toolkit//MATLAB/PC_Linux_GCC_MATLAB9.x_64bit/packages/mice.tar.Z
tar -xvf mice.tar.Z
cd mice && csh makeall.csh
```

### Get Kernels

`cd` to kernel directory.

```
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/sat441.bsp
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls.pc
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc
```