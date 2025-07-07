from dataclasses import dataclass
import numpy as np
import xarray as xr
from dataclasses import dataclass
from datetime import datetime
from zoneinfo import ZoneInfo
import rocketpy
from typing import Union
import pandas

import models

G = 9.80665

def drogue_trigger(p, h, y):
    # Check if rocket is going down, i.e. if it has passed the apogee
    vertical_velocity = y[5]
    # Return true to activate parachute once the vertical velocity is negative
    return vertical_velocity < 0 # False if you want ballistic

def main_trigger(p,h,y):
    vertical_velocity = y[5]
    return y[5] < 0 and h <= 400 # False if you want ballistic

def filter_flight_data(flight: rocketpy.Flight) -> models.Flight:
    return models.Flight(
        x_impact=flight.x_impact,
        y_impact=flight.y_impact,
    )



def construct_env(
    r: models.RocketParams,          # carries launch_time, elevation, etc.
    w: models.Weather,               # surface observation (speed, dir, T, p)
    time: datetime,
    lat: float,
    lon: float,
    climatology_file: str = "inputs/MC_env.nc",
) -> rocketpy.Environment:
    """
    Build a *single* RocketPy Environment for the rocket described in `r`
    and the surface weather `w`, by blending the MC_env.nc climatology
    with the current MET observation.
    """
    ds = xr.open_dataset(climatology_file)
    ts = np.datetime64(time.astimezone(ZoneInfo("UTC")))
    clim = ds.sel(time=ts, latitude=lat, longitude=lon, method="nearest")

    height         = (clim["z"].values / G).astype(float)       # m
    pressure_pa    = (clim["level"].values * 100.0).astype(float)  # Pa
    T_profile_K    = clim["t"].values                           # K

    T_obs_K        = w.temperature + 273.15                     # C → K
    T_surf_clim_K  = np.interp(0.0, height, T_profile_K)
    dT             = T_obs_K - T_surf_clim_K
    temperature_profile = [
        (float(h), float(T + dT)) for h, T in zip(height, T_profile_K)
    ]

    θ     = np.deg2rad(w.wind_direction)        # deg → rad (from North)
    u0    = -w.wind_speed * np.sin(θ)            # meteorological → u, v
    v0    = -w.wind_speed * np.cos(θ)
    wind_u_profile = [(float(h), float(u0)) for h in height]
    wind_v_profile = [(float(h), float(v0)) for h in height]

    # pressure profile (no anomaly applied)
    pressure_profile = [(float(h), float(p)) for h, p in zip(height, pressure_pa)]

    env = rocketpy.Environment(
        max_expected_height=getattr(r, "max_expected_height", 12_000),
        latitude   = lat,
        longitude  = lon,
        elevation  = getattr(r, "launch_elevation", 20.0)  # m ASL
    )
    env.set_date(r.launch_time, timezone="Europe/Oslo")
    env.set_atmospheric_model(
        type        = "custom_atmosphere",
        pressure    = pressure_profile,
        temperature = temperature_profile,
        wind_u      = wind_u_profile,
        wind_v      = wind_v_profile,
    )
    return env
