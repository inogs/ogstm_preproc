#!/usr/bin/env python3
import os
import glob
import numpy as np
import xarray as xr
import xesmf as xe
from scipy.spatial import cKDTree
from argparse import ArgumentParser
import yaml

'''
Re-grids OGSTM-BFM restart files from a source rectilinear grid to a curvilinear grid. 
Both have the same z-level, but the destination grid has some more deep levels.
Each depth level is re-gridded independently.
Values are extended at depth from the deepest y, x point
'''

# usage:
# source /leonardo/home/userexternal/ggalli00/miniconda3/etc/profile.d/conda.sh
# conda activate /leonardo_work/OGS_test2528_0/ggalli/RIVIERADE/PYENVS/regrid_env
# python regrid_restart.py -y regrid_restart.yaml

def load_parameters():
    parser = ArgumentParser()
    parser.add_argument(
        "--yamlfile", "-y",
        required=True,
        help="path to file with all parameters"
    )
    args = parser.parse_args()
    yamlfile = args.yamlfile
    # yamlfile = 'analyse_fluxes.yaml'
    with open(yamlfile) as f:
        Params = yaml.safe_load(f)
    return Params

Params = load_parameters()

# ============================================================
# USER SETTINGS
# ============================================================

# Directory containing source OGSTM-BFM restart files, files pattern, coordinates
#SOURCE_DIR = "/leonardo_scratch/fast/RIVRD_OGS/V13C_QUID/IC/"
SOURCE_DIR = Params['Input_files']['path']
FILE_PATTERN = Params['Input_files']['pattern']
SRC_LON = Params['Input_files']['lon_var']
SRC_LAT = Params['Input_files']['lat_var']
SRC_DEPTH = Params['Input_files']['depth_var']

# Directory where regridded restart files will be written
#OUTPUT_DIR = "/leonardo_scratch/fast/RIVRD_OGS/RST/"
OUTPUT_DIR = Params['Output_dir']

# NEMO destination mesh file and T-grid mask variable, destination grid coordinate names
#DEST_MESH_FILE = "/leonardo_work/RIVRD_OGS/OGSTM-BFM-SETUP/MASKS/med24v1_domain_cfg.nc"
#DEST_MASK_VAR = "tmask"
#DEST_LON = "glamt"
#DEST_LAT = "gphit"
DEST_MESH_FILE = Params['Destination_mesh']['file']
DEST_MASK_VAR = Params['Destination_mesh']['mask_var']
DEST_LON = Params['Destination_mesh']['lon_var']
DEST_LAT = Params['Destination_mesh']['lat_var']

# Source mask file and T-grid mask variable
#SRC_MASK_FILE = "/leonardo_scratch/fast/RIVRD_OGS/V13C_QUID/MASKS/meshmask.nc"
#SRC_MASK_VAR = "tmask"
#SRC_LON = "nav_lon"
#SRC_LAT = "nav_lat"
#SRC_DEPTH = "nav_lev"
SRC_MASK_FILE = Params['Source_mesh']['file']
SRC_MASK_VAR = Params['Source_mesh']['mask_var']

# Regridding method
REGRID_METHOD = "nearest_s2d"

# Minimum valid interpolation coverage.
# 1.e-12 avoids division by zero.
MIN_COVERAGE = 1.e-12

# Whether to overwrite existing output files
OVERWRITE = True

# ============================================================
# CREATE OUTPUT DIRECTORY
# ============================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================
# READ DESTINATION GRID
# ============================================================

print("Reading destination NEMO grid...")

with xr.open_dataset(DEST_MESH_FILE) as ds_mesh:

    dst_lon = ds_mesh[DEST_LON].squeeze().load()
    dst_lat = ds_mesh[DEST_LAT].squeeze().load()

    # NEMO 3D T-grid mask: 1 = ocean, 0 = land
    dst_mask3d = ds_mesh[DEST_MASK_VAR].squeeze().load()

    dst_depth = ds_mesh["nav_lev"].load()

# Mask dimensions are (z, y, x)
dst_mask3d = dst_mask3d.transpose(*dst_mask3d.dims[-3:])
dst_wet3d = dst_mask3d.values > 0.5
dst_mask = dst_mask3d.isel({dst_mask3d.dims[0]: 0})

dst_grid = xr.Dataset({
    "lon": (("y", "x"), dst_lon.values),
    "lat": (("y", "x"), dst_lat.values),
    "mask": (("y", "x"), dst_mask.values.astype(np.int32)),
})

NZ_DST = dst_wet3d.shape[0]

print("Destination grid shape:", dst_grid["lon"].shape, "levels:", NZ_DST)


# ============================================================
# READ SOURCE GRID
# ============================================================

restart_files = sorted(glob.glob(os.path.join(SOURCE_DIR, FILE_PATTERN)))

if not restart_files:
    raise FileNotFoundError("No source restart files found!")

print("Number of restart files:", len(restart_files))

print("Reading source grid from first restart...")

with xr.open_dataset(restart_files[0]) as ds_src:

    src_lon_2d = ds_src[SRC_LON].load()
    src_lat_2d = ds_src[SRC_LAT].load()

    src_depth = ds_src[SRC_DEPTH].load()

    # Extract 1D coordinates from the rectilinear grid
    # nav_lon is assumed constant along y
    # nav_lat is assumed constant along x

    src_lon = src_lon_2d.isel(y=0)
    src_lat = src_lat_2d.isel(x=0)

    # Identify the first biogeochemical tracer
    tracer_candidates = [v for v in ds_src.data_vars if v not in [SRC_LON, SRC_LAT, SRC_DEPTH]]

    if not tracer_candidates:
        raise ValueError("No tracer variables found!")

    first_tracer = ds_src[tracer_candidates[0]]


# read source mask
with xr.open_dataset(SRC_MASK_FILE) as ds_mask:

    src_wet = ds_mask[SRC_MASK_VAR].squeeze().load()

src_wet = src_wet.transpose(*src_wet.dims[-3:]).values > 0.5

NZ_SRC = src_wet.shape[0]
NY_SRC, NX_SRC = src_wet.shape[1:]

if NZ_SRC > NZ_DST:
    raise ValueError("Source grid has more levels than destination grid!")

# Both grids must share the same z-levels over the source depth range
if not np.allclose(src_depth.values, dst_depth.values[:NZ_SRC], rtol=1.e-4):
    raise ValueError("Source and destination depth levels differ over the common range!")

# Construct source grid. The regridder mask is the surface mask: any point
# wet at depth is also wet at the surface. Level-dependent land is handled
# in regrid_tracer.
src_grid_2d = xr.Dataset({
    "lon": ("x", src_lon.values),
    "lat": ("y", src_lat.values),
    "mask": (("y", "x"), src_wet[0].astype(np.int32)),
})

print("Source grid shape:", src_wet.shape)


# ============================================================
# CREATE XESMF REGRIDDER
# ============================================================

print("Creating xESMF regridder...")

regridder = xe.Regridder(
    src_grid_2d,
    dst_grid,
    method=REGRID_METHOD,
    periodic=False,
    unmapped_to_nan=True,
    filename="weights_bfm.nc",
    reuse_weights=False)

print("Regridder successfully created.")


# ============================================================
# REGRID FUNCTION
# ============================================================

def to_xyz(lon, lat):
    """Unit-sphere cartesian coordinates, so that KD-tree distances follow great circles."""
    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)
    return np.column_stack([
        (np.cos(lat) * np.cos(lon)).ravel(),
        (np.cos(lat) * np.sin(lon)).ravel(),
        np.sin(lat).ravel(),
    ])


src_lon2d, src_lat2d = np.meshgrid(src_lon.values, src_lat.values)
dst_xyz = to_xyz(dst_grid["lon"].values, dst_grid["lat"].values)

# KD-trees of the wet source points of each level (built lazily, reused)
_src_trees = {}


def src_level_tree(k):
    if k not in _src_trees:
        iy, ix = np.nonzero(src_wet[k])
        if iy.size == 0:
            _src_trees[k] = None
        else:
            tree = cKDTree(to_xyz(src_lon2d[iy, ix], src_lat2d[iy, ix]))
            _src_trees[k] = (tree, iy, ix)
    return _src_trees[k]


def fill_from_nearest_source(field, data, k, gap):
    """
    Fill the destination points `gap` (2D boolean) of level k with the value of
    the nearest wet source point of the SAME level.
    """
    entry = src_level_tree(k)

    if entry is None:
        return

    tree, iy, ix = entry
    _, nn = tree.query(dst_xyz[gap.ravel()])
    field[gap] = data[k][iy[nn], ix[nn]]


def fill_from_nearest_dest(field2d, wet2d):
    """
    Fill remaining NaNs at wet points of a 2D destination level with the nearest
    valid destination value of the same level.
    """
    valid = wet2d & np.isfinite(field2d)
    gap = wet2d & ~np.isfinite(field2d)

    if not gap.any() or not valid.any():
        return

    tree = cKDTree(dst_xyz[valid.ravel()])
    _, nn = tree.query(dst_xyz[gap.ravel()])
    field2d[gap] = field2d[valid][nn]


def regrid_3d(data):
    """
    Regrid one 3D field (z, y, x) defined on the source grid to the full
    destination grid (z_dst, y, x).

    Every level is regridded independently with horizontal weights only. Source
    land (NaN, 0 or fill values) is excluded through the source mask, so it
    cannot contaminate destination wet points. Destination wet points that
    cannot be reached by valid source values of their level are filled from the
    nearest wet source point of that level. Destination levels below the source
    domain are filled with the deepest valid value of the same (y, x) column.
    """

    valid = src_wet & np.isfinite(data)
    data_clean = np.where(valid, data, 0.0)

    # Weighted regridding, normalised by the weights of valid points only
    numerator = regridder(data_clean)
    denominator = regridder(valid.astype(np.float64))

    with np.errstate(invalid="ignore", divide="ignore"):
        result = np.where(denominator > MIN_COVERAGE, numerator / denominator, np.nan)

    out = np.full((NZ_DST,) + result.shape[1:], np.nan, dtype=np.float64)
    out[:NZ_SRC] = result
    out[~dst_wet3d] = np.nan

    # Same-level fill from the nearest wet source point
    for k in range(NZ_SRC):
        gap = dst_wet3d[k] & ~np.isfinite(out[k])
        if gap.any():
            fill_from_nearest_source(out[k], np.where(valid, data, np.nan), k, gap)

    # Vertical extrapolation of the deepest valid value (also used for any
    # column that has no valid value at some level)
    finite = np.isfinite(out)
    last = np.maximum.accumulate(np.where(finite, np.arange(NZ_DST)[:, None, None], -1), axis=0)
    filled = np.take_along_axis(out, np.maximum(last, 0), axis=0)
    filled[last < 0] = np.nan
    out = np.where(dst_wet3d, filled, np.nan)

    # Last resort: nearest valid value on the same destination level
    for k in range(NZ_DST):
        fill_from_nearest_dest(out[k], dst_wet3d[k])

    if np.isnan(out[dst_wet3d]).any():
        raise RuntimeError("Some destination wet points could not be filled!")

    return out


def regrid_tracer(data):
    """
    Regrid a biogeochemical tracer.

    Parameters
    ----------
    data : xarray.DataArray
        Tracer with dimensions including (z, y, x) and optionally leading
        dimensions such as time.

    Returns
    -------
    xarray.DataArray
        Regridded tracer with NZ_DST levels, NaN on destination land.
    """

    lead = [d for d in data.dims if d not in ("z", "y", "x")]
    data = data.transpose(*lead, "z", "y", "x")
    values = data.values

    if values.shape[-3:] != src_wet.shape:
        raise ValueError(f"Unexpected tracer shape {values.shape[-3:]}, expected {src_wet.shape}")

    out = np.empty(values.shape[:-3] + (NZ_DST,) + dst_grid["lon"].shape, dtype=np.float64)

    for idx in np.ndindex(*values.shape[:-3]):
        out[idx] = regrid_3d(values[idx])

    result = xr.DataArray(out.astype(data.dtype), dims=data.dims, attrs=data.attrs.copy())

    return result


# ============================================================
# PROCESS RESTART FILES
# ============================================================

for filename in restart_files:

    basename = os.path.basename(filename)

    output_file = os.path.join(OUTPUT_DIR, basename)

    if os.path.exists(output_file) and not OVERWRITE:
        print("Skipping existing file:", basename)
        continue

    print("\nProcessing:", basename)

    with xr.open_dataset(filename) as ds:

        # Skip files that are not on the source grid (e.g. already regridded)
        if ds.sizes.get("y") != NY_SRC or ds.sizes.get("x") != NX_SRC:
            print(f"  WARNING: skipping {basename}: horizontal shape "
                  f"({ds.sizes.get('y')}, {ds.sizes.get('x')}) differs from the "
                  f"source grid ({NY_SRC}, {NX_SRC}).")
            continue

        # Load source dataset into memory
        ds = ds.load()

        # Create output dataset
        ds_out = xr.Dataset(attrs=ds.attrs.copy())

        # ----------------------------------------------------
        # Copy time and vertical coordinates
        # ----------------------------------------------------

        if "time" in ds.coords:
            ds_out = ds_out.assign_coords(time=ds["time"])

        # Destination vertical levels (includes the additional deep levels)
        ds_out["nav_lev"] = ("z", dst_depth.values.astype(ds["nav_lev"].dtype))
        ds_out["nav_lev"].attrs = ds["nav_lev"].attrs.copy()

        # ----------------------------------------------------
        # Add destination horizontal coordinates
        # ----------------------------------------------------

        ds_out["nav_lon"] = dst_grid["lon"]
        ds_out["nav_lat"] = dst_grid["lat"]

        ds_out["nav_lon"].attrs = ds[SRC_LON].attrs.copy()
        ds_out["nav_lat"].attrs = ds[SRC_LAT].attrs.copy()

        # ----------------------------------------------------
        # Regrid all biogeochemical tracers
        # ----------------------------------------------------

        for varname in ds.data_vars:

            if varname in [SRC_LON, SRC_LAT, SRC_DEPTH]:
                continue

            variable = ds[varname]

            # Only regrid variables with horizontal dimensions
            if not {"z", "y", "x"}.issubset(variable.dims):
                continue

            print("  Regridding:", varname)

            # Regrid tracer
            ds_out[varname] = regrid_tracer(variable)

        # ----------------------------------------------------
        # Write output
        # ----------------------------------------------------

        print("  Writing:", output_file)

        encoding = {}

        for varname in ds_out.data_vars:

            if varname in ["nav_lev", "nav_lon", "nav_lat"]:
                continue

            if np.issubdtype(ds_out[varname].dtype, np.floating):

                encoding[varname] = {
                    "_FillValue": ds_out[varname].dtype.type(1.e20),
                    "zlib": True,
                    "complevel": 4
                }

        ds_out.to_netcdf(output_file, encoding=encoding)

    print("  Done.")

# ============================================================
# CLEANUP
# ============================================================

print("\nAll restart files successfully regridded.")

