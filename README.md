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


## Running code.

Call `lofi_search` followed by any number of input files.

e.g. 

```matlab
lofi_search    TEE_basics    hpc_conf
```

Where `TEE_basics` contains mission info, and `hpc_conf` contains site specific info.


### Note

site specific setup should be moved outside of config eventually.

## Using Git

### Checking status 

```matlab
!git status
```

### Before you make any changes to the code

```matlab
!git pull
```

There shouldn't ever be issues with this command. Seek help if there are.

### When you are done making changes

```matlab
!git add <files you want recorded>
!git commit -m "description of changes"
!git push
```

If there are no files in your working directory that you don't want recorded (and arn't listed in `.gitignore`) you can use `!git add --all` instead of listing files.

## Kernels

Note, paths in meta kernel files are relative to current working directory, not kernel file. 
Therefore, scripts should always be run from same directory else meta-kernels will not work.