import numpy as np
import pandas as pd # type: ignore
from dataclasses import asdict
from pathlib import Path
from typing import Any


import models


xlsx = Path("input_values.xlsx")
if not xlsx.exists():
    raise FileNotFoundError(f"{xlsx} not found")

df = pd.read_excel(xlsx, index_col=1)
df.columns = df.iloc[0]  # promote first row → column names
df = df.iloc[1:]  # drop the old header row


def v(row: str, pos: int = 1) -> Any:
    """Shortcut: value in the *second* column of the given row."""
    return df.loc[row].iloc[pos]


# ───────────────────────────────── baseline analysis settings
ANALYSIS_PARAMETERS: dict[str, Any] = {
    # Mass details
    "rocket_length": (v("rocket_length"), 0.05),
    "radius": (v("rocket_radius"), 0.001),
    "cg": (v("center_gravity"), 0.05),
    "dry_mass": (v("rocket_mass"), 1),
    "Ixx": (v("inertia_xx"), 0.05),
    "Iyy": (v("inertia_yy"), 3),
    "Izz": (v("inertia_zz"), 3),
    "Ixy": (v("inertia_xy"), 0.1),
    "Ixz": (v("inertia_xz"), 0.05),
    "Iyz": (v("inertia_yz"), 0.05),
    # Engine and tank details
    "burnout_time": (9, 1),
    "tank_length": (v("fuel_length"), 0),
    "fuel_tank_radius": (v("fuel_radius"), 0),
    "ox_tank_length": (v("ox_length"), 0),
    "ox_tank_radius": (v("ox_eff_radius"), 0),
    "n2_tank_radius": (v("n2_radius"), 0),
    "n2_tank_length": (v("n2_length"), 0),
    "fuel_tank_position": (v("fuel_position"), 0.05),
    "ox_tank_position": (v("ox_position"), 0.05),
    "n2_tank_position": (v("n2_position"), 0.05),
    "fuel_mass": (v("fuel_mass"), 0),
    "n2_volume": (v("n2_volume"), 0),
    "ox_mass": (v("ox_mass"), 0),
    "of": (v("OF_ratio"), 0.2),
    "fuel_pressure": (v("fuel_pressure"), 1.0e6),
    "ox_pressure": (v("ox_pressure"), 1.0e6),
    "n2_pressure": (v("n2_pressure"), 1.0e6),
    "ox_temp": (v("ox_temp"), 0),
    "fuel_temp": (v("fuel_temp"), 0),
    "ambient_temp": (v("ambient_temp"), 0),
    "ethanol_perc": (v("ethanol_perc"), 0),
    "water_perc": (v("water_perc"), 0),
    "mdot": (v("Massflowrate"), 0.1),
    # Aerodynamic details
    "nozzle_position": (0, 0.01),
    "power_off_drag": (1, 0),
    "power_on_drag": (1, 0),
    "nose_length": (v("nose_length"), 0.001),
    "fin_span": (v("fin_span"), 0.0005),
    "fin_root_chord": (v("rootchord"), 0.0005),
    "fin_tip_chord": (v("tipchord"), 0.0005),
    "beta": (v("fin_beta"), 1),
    "sweep_length": (v("fin_span") / np.tan(np.deg2rad(v("fin_beta"))), 0.001),
    "fin_position": (v("fin_position"), 0.05),
    "inclination": (v("inclination"), 1),
    "heading": (v("heading"), 2),
    "rail_length": (v("rail_length"), 0.005),
    "ensemble_member": list(range(10)),
    "drogue_Cd_s": (v("drogue_cds"), 0.01 * v("drogue_cds")),
    "drogue_lag": (v("drogue_total_lag"), 0.01 * v("drogue_total_lag")),
    "main_Cd_s": (v("main_cds"), 0.01 * v("main_cds")),
    "main_lag": (v("main_total_lag"), 0.01 * v("main_total_lag")),
}


def _replace_mean(old: tuple, new_val: float) -> tuple:
    """
    Keep the existing stdev but replace the mean with `new_val`.
    """
    return (new_val, *old[1:])


def get_parameters(r: models.RocketParams) -> dict[str, Any]:
    """
    Return a fresh dict equal to ANALYSIS_PARAMETERS except that any
    attribute present on `r` overrides the *mean* value.
    """
    params = ANALYSIS_PARAMETERS.copy()

    r_values = asdict(r)

    for key, val in r_values.items():
        if key in params and isinstance(params[key], tuple):
            params[key] = _replace_mean(params[key], val)

    return params
