import mpi4py
import blend_ens_align
import xarray as xr
import numpy as np
import datetime
import yaml
import sys

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

s = datetime.datetime.now()
if rank == 0:
    print("*"*40)
    print(f" BEGIN test_ens_align - {s.strftime('%Y-%m-%d %H:%M:%S')}")
    print("*"*40)

select_models = sys.argv[1]
select_models = select_models.split()
out_zarr = sys.argv[2]

"""
In operation, the QMD model forecast files will be less clean, so we will rely on 
this function from precip prodgen to grab only the QMD forecasts
"""
def collect_augment_das(ds: xr.Dataset, model_names: []) -> list:

    # Use center grid point to check if member is available.
    xchk = int(ds.sizes['xa'] / 2)
    ychk = int(ds.sizes['ya'] / 2)

    # Iterate over model DataArrays
    da_models = []
    for k, v in ds.data_vars.items():
        if "fcst" in k:

            if "qmd" in k:
                da = v.squeeze(dim=['time', 'lead_time', 'nstencil'])
                if (len(model_names) > 0) and (da.name not in model_names):
                    continue
                shortname = da.name.split("_")[0]
                kind = 'QMD'
            else:
                continue

            da = da.rename({da.dims[0]: 'nmembers'})
            da = da.assign_coords({'nmembers': [k.split('_')[0] for i in range(da.nmembers.size)]})
            da.attrs['nmembers_expected'] = da.nmembers.size

            # Check for model availability here...
            miss_count = np.count_nonzero(np.isnan(da[:, ychk, xchk].values))
            if miss_count > 0:
                if miss_count == da.nmembers.size:
                    print(f" ****Model {da.name} is completely missing. Will not be included in list of DataArrays.")
                    continue # Next model up...
                elif miss_count < da.nmembers.size:
                    print(f" ****Model {da.name} is missing {miss_count} members. Will trim the DataArray.")
                    valid_mask = ~np.isnan(da[:, ychk, xchk].values)
                    da = da.isel(nmembers=valid_mask)

            da.attrs['nmembers_available'] = da.nmembers.size
            da = da.astype("float32")
            da_models.append(da)

    return da_models

ens_data = xr.open_dataset('test_data/blend.t00z.model_qmd_fcst.precip06.f006.co.2p5.nc')

lat = ens_data.latitude.data
lon = ens_data.longitude.data
ny, nx = np.shape(lat) # --> remember that Fortran index = Python+1 AND x,y indices flipped
ifhr = ens_data.lead_time.values[0]
dx, dy = 2500., 2500.
    
if rank == 0:
    print(f"[PYTHON]: Loading data...", flush=True)
    print(f"           Using models: {select_models}", flush=True)

ens_das = collect_augment_das(ens_data, select_models)

member_loc = dict()
for k in ens_das:
    member_loc[k.name] = k.nmembers.size

### fortran expects forecast data in shape of ensfcst(nx,ny,membknt)

ensfcst = np.zeros((nx, ny, sum(member_loc.values())))
ind = 0
for da in ens_das:
    nmem = da.nmembers.size
    # set any NaNs to missing value
    ensfcst[:,:,ind:ind+nmem] = np.nan_to_num(np.asfortranarray(da.data.transpose(2, 1, 0)), nan=-999.)
    ind += nmem

with open("configs/blend.qmd.precip24.ens_align.yaml", "r") as f:
    config = yaml.safe_load(f)


nbaksmth,nshfsmth,applyshft = config['shift_const']['nbaksmth'],config['shift_const']['nshfsmth'],config['shift_const']['applyshft']
slnratio0h,slnratio48h = config['shift_const']['slnratio0h'], config['shift_const']['slnratio48h']
noutsmth,minkdat,minkdratio = config['shift_const']['noutsmth'],config['shift_const']['minkdat'],config['shift_const']['minkdratio']

hrzlap = config['shift_zone']['hrzlap']
iborder = config['shift_zone']['iborder']
izsize = config['shift_zone']['izsize']
jborder = config['shift_zone']['jborder']
jzsize = config['shift_zone']['jzsize']
loopstep,procspg = config['shift_zone']['loopstep'], config['shift_zone']['procspg']

wgtvar,thresh_flag,threshvar = config['shift_wg']['wgtvar'], config['shift_wg']['thresh_flag'],config['shift_wg']['threshvar']


gauss_sigma,patch_nx,patch_ny,ovx,ovy,filt_min = (config['lpm_const']['gauss_sigma'], config['lpm_const']['patch_nx'],
                                                  config['lpm_const']['patch_ny'], config['lpm_const']['ovx'],
                                                  config['lpm_const']['ovy'], config['lpm_const']['filt_min'])

if rank == 0:
    print(f"[PYTHON]: Sending to align_blend_ens.ensalign", flush=True)
    print(f"           N members: {ind}", flush=True)


shift_x, shift_y, ensfcst_shift,ensfcst_shift_pm, ensfcst_shift_lpm = blend_ens_align.ensalign(ensfcst,dy,dx,ifhr,
                                   nbaksmth,nshfsmth,applyshft,slnratio0h,slnratio48h,noutsmth,minkdat,minkdratio,
                                   hrzlap,iborder,izsize,jborder,jzsize,loopstep,procspg,
                                   wgtvar,thresh_flag,threshvar,
                                   patch_nx,patch_ny,ovx,ovy,gauss_sigma,filt_min)

if rank == 0:

    results = xr.Dataset(
                data_vars=dict(
                    x_shift_field=(["shift_pass", "member", "y", "x"], np.ascontiguousarray(shift_x.transpose(3,2, 1, 0))),
                    y_shift_field=(["shift_pass", "member", "y", "x"], np.ascontiguousarray(shift_y.transpose(3,2, 1, 0))),
                    qpf_aligned=(["shift_pass", "member", "y", "x"], np.ascontiguousarray(ensfcst_shift.transpose(3,2, 1, 0))),
                    qpf_aligned_pm=(["shift_pass","y", "x"], np.ascontiguousarray(ensfcst_shift_pm.transpose(2, 1, 0))),
                    qpf_aligned_lpm=(["shift_pass","y", "x"], np.ascontiguousarray(ensfcst_shift_lpm.transpose(2, 1, 0)))
                ),
                coords=dict(
                    longitude=(['y', 'x'], lon),
                    latitude=(['y', 'x'], lat),
                    time=ens_data.time.data,
                    lead_time=ens_data.lead_time.data,
            
                ),
                attrs=dict(description="Ensemble alignment results"),
            )

    print(f"[PYTHON]: Saving results to {out_zarr}")
    results.to_zarr(out_zarr, mode='w')

    
    print("*"*40)
    print(f" FINISHED in {(datetime.datetime.now() - s).total_seconds():.2f} seconds")
    print("*"*40)