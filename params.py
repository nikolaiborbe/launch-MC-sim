import numpy as np
import pandas as pd

df = pd.read_excel("input_values.xlsx", index_col=1)

# --- tidy the columns exactly as before ---
headers = df.iloc[0].values          # old first row becomes column names
df.columns = headers
df = df.iloc[1:]                     # drop the header row itself, keep data

# Helper: second *column by position* of any row
def v(row_name, pos=1):
    return df.loc[row_name].iloc[pos]

ANALYSIS_PARAMETERS = {
    # Mass Details
    "rocket_length": (v("rocket_length"), 0.05),
    "radius":        (v("rocket_radius"), 0.001),
    "cg":            (v("center_gravity"), 0.05),
    "dry_mass":      (v("rocket_mass"), 1),
    "Ixx":           (v("inertia_xx"), 0.05),
    "Iyy":           (v("inertia_yy"), 3),
    "Izz":           (v("inertia_zz"), 3),
    "Ixy":           (v("inertia_xy"), 0.1),
    "Ixz":           (v("inertia_xz"), 0.05),
    "Iyz":           (v("inertia_yz"), 0.05),

    # ENGINE AND TANK DETAILS
    "burnout_time":          (9, 1),
    "tank_length":           (v("fuel_length"), 0),
    "fuel_tank_radius":      (v("fuel_radius"), 0),
    "ox_tank_length":        (v("ox_length"), 0),
    "ox_tank_radius":        (v("ox_eff_radius"), 0),
    "n2_tank_radius":        (v("n2_radius"), 0),
    "n2_tank_length":        (v("n2_length"), 0),
    "fuel_tank_position":    (v("fuel_position"), 0.05),
    "ox_tank_position":      (v("ox_position"), 0.05),
    "n2_tank_position":      (v("n2_position"), 0.05),
    "fuel_mass":             (v("fuel_mass"), 0),
    "n2_volume":             (v("n2_volume"), 0),
    "ox_mass":               (v("ox_mass"), 0),
    "of":                    (v("OF_ratio"), 0.2),
    "fuel_pressure":         (v("fuel_pressure"), 10e5),
    "ox_pressure":           (v("ox_pressure"), 10e5),
    "n2_pressure":           (v("n2_pressure"), 10e5),
    "ox_temp":               (v("ox_temp"), 0),
    "fuel_temp":             (v("fuel_temp"), 0),
    "ambient_temp":          (v("ambient_temp"), 0),
    "ethanol_perc":          (v("ethanol_perc"), 0),
    "water_perc":            (v("water_perc"), 0),
    "mdot":                  (v("Massflowrate"), 0.1),

    # Aerodynamic Details
    "nozzle_position": (0, 0.01),
    "power_off_drag":  (1, 0),
    "power_on_drag":   (1, 0),
    "nose_length":     (v("nose_length"), 0.001),
    "fin_span":        (v("fin_span"), 0.0005),
    "fin_root_chord":  (v("rootchord"), 0.0005),
    "fin_tip_chord":   (v("tipchord"), 0.0005),
    "beta":            (v("fin_beta"), 1),
    "sweep_length":    (v("fin_span") / np.tan(np.deg2rad(v("fin_beta"))), 0.001),
    "fin_position":    (v("fin_position"), 0.05),
    "inclination":     (v("inclination"), 1),
    "heading":         (v("heading"), 2),
    "rail_length":     (v("rail_length"), 0.005),
    "ensemble_member": list(range(10)),
    "drogue_Cd_s":     (v("drogue_cds"), 0.01 * v("drogue_cds")),
    "drogue_lag":      (v("drogue_total_lag"), 0.01 * v("drogue_total_lag")),
    "main_Cd_s":       (v("main_cds"), 0.01 * v("main_cds")),
    "main_lag":        (v("main_total_lag"), 0.01 * v("main_total_lag")),
}
