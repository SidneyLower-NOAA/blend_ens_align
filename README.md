# Spatially Aligned Mean for NBM

Adapted from source at: https://github.com/hunter3789/ensemble-align/tree/main

## Updates from source
`f12py` is used to wrap the source Fortran subroutines and modules to build `blend_ens_align` containing only the blocks relevant to alignment and transfers IO responsibility to Python via `xarray`. Alignment follows `alignopt=2` in the source, meaning ensemble members are aligned to ensemble probability matched mean (PM-mean) location. 

The `src/` directory contains the consolidated subroutines and modules `ensalign` depends on. 

- align_mod: common variables module, adapted from align.inc
- field_ops: grid editing methods and other utlities
- gs_smooth: gaussian smoothing routines
- mpi_utils: MPI broadcast adaptors for large arrays - Fortran's answer to xarray chunking and lazy loading
- mrgrnk: merge-sort ranking module
- pmm_lpm: Global and localized probability matched mean routines
- shift_ops: shift/alignment vector routines
- ens_align_blend: primary subroutine adapted from `ensalign` main program

Variables previously set by namelists are now included in yaml files in `configs/`

## Build
To build the Fortran modules with `f2py`, follow the guide in `build_ens_align.sh` -- the exact compiler suite listed there is specific to Ursa so use whatever stack you prefer. 

## How to run
`blend_ens_align.ensalign` requires the following inputs:

```
ensfcst : input rank-3 array('f') with bounds (py_nx,py_ny,nmembers)
py_dy : input int
py_dx : input int
ifhr : input int
py_nbaksmth : input int
py_nshfsmth : input int
py_applyshft : input int
py_slnratio0h : input float
py_slnratio48h : input float
py_noutsmth : input int
py_minkdat : input float
py_minkdratio : input float
py_hrzlap : input float
py_iborder : input rank-1 array('i') with bounds (py_nshfpass)
py_izsize : input rank-1 array('i') with bounds (py_nshfpass)
py_jborder : input rank-1 array('i') with bounds (py_nshfpass)
py_jzsize : input rank-1 array('i') with bounds (py_nshfpass)
py_loopstep : input rank-1 array('i') with bounds (py_nshfpass)
py_procspg : input rank-1 array('i') with bounds (py_nshfpass)
py_wgtvar : input float
py_thresh_flag : input int
py_threshvar : input float
py_patch_nx : input int
py_patch_ny : input int
py_ovx : input int
py_ovy : input int
py_gauss_sigma : input int
py_filt_min : input int
```
And returns
```
xshiftmn : rank-4 array('f') with bounds (py_nx,py_ny,nmembers,py_nshfpass)
yshiftmn : rank-4 array('f') with bounds (py_nx,py_ny,nmembers,py_nshfpass)
ensfcst_shf : rank-4 array('f') with bounds (py_nx,py_ny,nmembers,py_nshfpass)
ensshfpm : rank-3 array('f') with bounds (py_nx,py_ny,py_nshfpass)
ensshflpm : rank-3 array('f') with bounds (py_nx,py_ny,py_nshfpass)
```

To run as an interactive job on Ursa with a sub-set of the full NBM suite:

```
export models="gefs_qmdfcst hrrrco_qmdfcst refs001_qmdfcst refs002_qmdfcst refs003_qmdfcst refs004_qmdfcst refs005_qmdfcst"
export outfile=ens_align_results.zarr

mpirun -np $NPROCS python test_ens_align.py 
```

To submit to SLURM scheduler:
```
sbatch run_ens_align.job
```

`test_ens_align.py` reads in an NBM QPF06 NetCDF file and extracts only the requested quantile mapped model arrays ('X_qmdfcst'). Numpy C-ordered arrays are transposed to align with Fortran style conventions. MPI initialization should occur in the runtime python script. Once `blend_ens_align.ensalign` finishes, xarray is used to save the output in `zarr` format.
