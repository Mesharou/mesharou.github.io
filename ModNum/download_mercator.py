#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download monthly subsets of GLORYS12V1 ocean reanalysis from Copernicus Marine Service
https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description

Requires Copernicus Marine Python API:
https://pypi.org/project/copernicusmarine/

Memory optimization:
- Avoid xarray reading + full rewrite just to correct time coordinates.
- Patch the `time` variable in-place using netCDF4 (edits only the coordinate).
"""

import os
import argparse
from datetime import datetime, timedelta
from dateutil import rrule, relativedelta
import copernicusmarine
from pathlib import Path

from Modules.grid_tools import define_subset_extent, find_croco_grd
from Modules.config_tools import load_ini, check_extent


def check_credentials(username_arg=None, password_arg=None):
    """Check if Copernicus Marine credentials exist or use provided login/password"""
    username = username_arg or os.getenv("CMEMS_USERNAME")
    password = password_arg or os.getenv("CMEMS_PASSWORD")

    # Try login/password arguments
    if username and password:
        try:
            status = copernicusmarine.login(
                username=username, password=password, force_overwrite=True
            )
            if status:
                print("✅ Logged in successfully using CLI arguments.")
                return
            else:
                raise SystemExit("❌ Login failed: invalid username or password.")
        except Exception as e:
            print(f"❌ Login with provided credentials failed: {e}")
    else:
        creds_file = Path.home() / ".copernicusmarine" / ".copernicusmarine-credentials"
        if creds_file.exists():
            print("🔑 Credentials detected → no login required.")
            return

    # Fall back to interactive login
    print("⚠️  No credentials found.")
    choice = input("Do you want to log in now? [y/N]: ").strip().lower()
    if choice == "y":
        copernicusmarine.login()
        print("✅ Logged in successfully.")
    else:
        raise SystemExit("❌ No credentials available, exiting.")


def fix_time_inplace_nc(nc_path: str, mode_tag: str) -> None:
    """
    Patch the 'time' variable in-place using netCDF4.

    mode_tag:
      - "P1D-m": daily means -> shift time by +12h
      - "P1M-m": monthly means -> set time to mid-month

    This edits only the coordinate variable and avoids rewriting the full file.
    """
    try:
        from netCDF4 import Dataset, num2date, date2num
    except Exception as e:
        raise RuntimeError(
            "netCDF4 is required for in-place time fixing. "
            "Install it with: mamba install -c conda-forge netcdf4"
        ) from e

    with Dataset(nc_path, "r+") as nc:
        if "time" not in nc.variables:
            return

        tvar = nc.variables["time"]

        # netCDF CF time metadata
        units = getattr(tvar, "units", None)
        if not units:
            raise RuntimeError(f"'time' variable has no 'units' attribute in {nc_path}")
        calendar = getattr(tvar, "calendar", "standard")

        # Read only the time coordinate (small)
        tvals = tvar[:]
        dts = num2date(tvals, units=units, calendar=calendar)

        if mode_tag == "P1D-m":
            new_dts = [dt + timedelta(hours=12) for dt in dts]

        elif mode_tag == "P1M-m":
            new_dts = []
            for dt in dts:
                start = datetime(dt.year, dt.month, 1)
                next_month = start + relativedelta.relativedelta(months=1)
                mid = start + (next_month - start) / 2
                new_dts.append(mid)

        else:
            return

        # Write back only the time values
        tvar[:] = date2num(new_dts, units=units, calendar=calendar)
        nc.sync()


def main():
    parser = argparse.ArgumentParser(
        description="""
            Download monthly subsets of GLORYS12V1 from Copernicus Marine Service
            Allows geographic area, variable and depth selection
        """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "filePath",
        nargs="?",
        default="download_mercator.ini",
        type=str,
        help="""
Path to the config file (default: download_mercator.ini)

Example : 
    [Times]
    Ystart = 2013 
    Mstart = 1 
    Yend = 2013 
    Mend = 1 

    [Download_Options]
    use_grd_extent = True
    grd_extent_padding = 1.0
    custom_extent = [latmax, lonmin, latmin, lonmax]

    [Croco_Files]
    croco_files_dir = ../results/CROCO_FILES
    croco_grd_prefix = croco_grd

    [IBC_Input_Files]
    ibc_dir = ../data/DATA_GLORYS_Benguela_LR/
    ibc_prefix = GLORYS 

    [Mercator_Download]
    dataset = cmems_mod_glo_phy_my_0.083deg_P1D-m
    variables = ["thetao","so","uo","vo","zos"]
    depths = [0,50]
        """,
    )
    parser.add_argument(
        "--username",
        type=str,
        help="""
Copernicus Marine username. 
If this argument is not given, CMEMS_USERNAME environment variable is used.
If none work, credential file 
( ~/.copernicusmarine/.copernicusmarine-credentials )
will be checked.
If credential file is not found, fall back to interactive login
        """,
    )
    parser.add_argument(
        "--password",
        type=str,
        help="""
Copernicus Marine password. 
If this argument is not given, CMEMS_PASSWORD environment variable is used.
If none work, credential file 
( ~/.copernicusmarine/.copernicusmarine-credentials )
will be checked.
If credential file is not found, fall back to interactive login
        """,
    )

    args = parser.parse_args()
    config = load_ini(args.filePath)

    check_credentials(username_arg=args.username, password_arg=args.password)

    # Load configuration
    Ystart = config["Times"]["Ystart"]
    Mstart = config["Times"]["Mstart"]
    Yend = config["Times"]["Yend"]
    Mend = config["Times"]["Mend"]

    dataset_id = config["Mercator_Download"]["dataset"]
    output_dir = config["IBC_Input_Files"]["ibc_dir"]
    output_name = config["IBC_Input_Files"]["ibc_prefix"]

    use_grd_extent = config["Download_Options"]["use_grd_extent"]

    variables = config["Mercator_Download"]["variables"]
    depth_cfg = config["Mercator_Download"]["depths"]

    date_start = datetime(Ystart, Mstart, 1)
    date_end = datetime(Yend, Mend, 1)

    # Normalize variables
    if isinstance(variables, str):
        variables = [v.strip() for v in variables.split(",")]
    elif isinstance(variables, list):
        variables = [v.strip() for v in variables]
    else:
        raise ValueError(f"Invalid variables config: {variables}")

    # Parse depth selection
    if isinstance(depth_cfg, str):
        if depth_cfg.lower() == "all":
            depth_min, depth_max = None, None
        else:
            try:
                values = [float(x.strip()) for x in depth_cfg.split(",")]
                if len(values) != 2:
                    raise ValueError
                depth_min, depth_max = values
            except Exception:
                raise ValueError(
                    f"Invalid depth string: {depth_cfg}. Must be 'all' or 'min,max'."
                )
    elif isinstance(depth_cfg, (list, tuple)):
        if len(depth_cfg) != 2:
            raise ValueError(
                f"Invalid depth list/tuple: {depth_cfg}. Must have 2 elements."
            )
        depth_min, depth_max = float(depth_cfg[0]), float(depth_cfg[1])
    else:
        raise ValueError(f"Depth config must be 'all', list, or tuple. Got {depth_cfg}")

    os.makedirs(output_dir, exist_ok=True)

    # Compute area
    if use_grd_extent:
        dl = config["Download_Options"]["grd_extent_padding"]
        croco_dir = config["Croco_Files"]["croco_files_dir"]
        croco_grd_prefix = config["Croco_Files"]["croco_grd_prefix"]
        croco_grd = find_croco_grd(croco_dir, croco_grd_prefix)
        area = define_subset_extent(croco_grd, dl)
        latmax, lonmin, latmin, lonmax = area
    else:
        custom_extent = config["Download_Options"]["custom_extent"]
        check_extent(custom_extent)
        latmax, lonmin, latmin, lonmax = custom_extent

    # Main loop over months
    for current_month in rrule.rrule(rrule.MONTHLY, dtstart=date_start, until=date_end):
        start_str = current_month.strftime("%Y-%m-%d %H:%M:%S")
        end_str = (
            current_month
            + relativedelta.relativedelta(months=1)
            - timedelta(seconds=1)
        ).strftime("%Y-%m-%d %H:%M:%S")

        print("Download from ", start_str, " to ", end_str)

        part_filename = f"{output_name}_Y{current_month.year}M{current_month.month:02d}.nc"
        part_filepath = os.path.join(output_dir, part_filename)

        # Download subset
        copernicusmarine.subset(
            dataset_id=dataset_id,
            minimum_longitude=lonmin,
            maximum_longitude=lonmax,
            minimum_latitude=latmin,
            maximum_latitude=latmax,
            start_datetime=start_str,
            end_datetime=end_str,
            variables=variables,
            minimum_depth=depth_min,
            maximum_depth=depth_max,
            output_filename=part_filename,
            output_directory=output_dir,
        )

        # Fix times in-place (no full rewrite)
        if dataset_id.endswith("P1D-m"):
            print(
                "Data are daily averages: shift times from 0h to 12h "
                "(in-place correction due to bad dates in Copernicus subset API)"
            )
            fix_time_inplace_nc(part_filepath, "P1D-m")

        elif dataset_id.endswith("P1M-m"):
            print(
                "Data are monthly averages: change time from 1st of month to mid-month "
                "(in-place correction due to bad dates in Copernicus subset API)"
            )
            fix_time_inplace_nc(part_filepath, "P1M-m")

    print("Download completed.")


if __name__ == "__main__":
    main()