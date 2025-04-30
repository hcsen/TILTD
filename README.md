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

## Using SLURM

The file `submit.sh` is a Slurm Script, a bash script with additional meta-parameters required for submission to the NeSI job scheduler.
Number of cpus, runtime and memory can all be set in the header of the script.

### Commands

|           |                       |                                                                          |
| --------- | --------------------- | ------------------------------------------------------------------------ |
| `sbatch`  | `sbatch submit.sl`    | Submits the Slurm script `submit.sl`                                      |
| `sacct`   | `sacct`               | Displays all the jobs run by you that day.                               |
|           | `sacct -S YYYY-MM-DD` | Displays all the jobs run by you since given date.                 |
| `scancel` | `scancel 123456789`   | Cancels job *123456789*                                                   |

Whatever output would normally be printed to stdout will instead be written to an output file, the name of the output file can be set in your Slurm script, but by default it will be `slurm-<jobid>.out`

## Kernels

Note, paths in meta kernel files are relative to current working directory, not kernel file. 
Therefore, scripts should always be run from same directory else meta-kernels will not work.

As kernels are not recorded in the repo (too big), I reccomend making a `Kernels` directory and putting all kernels in there.