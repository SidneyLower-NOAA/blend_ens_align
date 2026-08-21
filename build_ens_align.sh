# Ursa install
module purge
module load intel-oneapi-compilers/2026.1.0
module load intel-oneapi-mpi/2021.18.0
module load netcdf-fortran/4.6.1
module load netcdf-c/4.9.2
module load cmake/3.30.2

# activate conda env
conda activate ens_align

# generate f2py signature file
python -m numpy.f2py ens_align_blend.f90 -h only_ens_align.pyf -m blend_ens_align --overwrite-signature

# compile with f2py
FC=mpif90 CC=mpicc python -m numpy.f2py -c --dep mpi only_ens_align.pyf \
  ens_align_blend.f90 field_ops.f90 gs_smooth.f90 mpi_utils.f90 mrgrnk.f90 pmm_lpm.f90 align_mod.f90 shift_ops.f90
### if debugging, add this line as well
  --f90flags="-g -O0 -fcheck=all -fbacktrace"
