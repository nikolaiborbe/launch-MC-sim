from params import ANALYSIS_PARAMETERS, df
import time
from utils import drogue_trigger, main_trigger
import pandas as pd
import pyfluids
from rocketpy import (
    Environment, Rocket, Flight, LiquidMotor, MassFlowRateBasedTank,
    CylindricalTank, Fluid
    
)
import numpy
from typing import Union

analysis_parameters: pd.DataFrame = ANALYSIS_PARAMETERS


def flight_settings(analysis_parameters: Union[tuple, float], total_number: int):
    i = 0
    while i < total_number:
        # Generate a flight setting
        flight_setting = {}
        for parameter_key, parameter_value in analysis_parameters.items():
            if type(parameter_value) is tuple:
                flight_setting[parameter_key] = numpy.random.normal(*parameter_value)
            else:
                flight_setting[parameter_key] = numpy.random.choice(parameter_value)

        # Update counter
        i += 1

        # Yield a flight setting
        yield flight_setting


def get_MC_sim_result(num_sims: int = 10) -> list[Flight]:
    MC_sim_result: list[Flight] = []

    for setting in flight_settings(analysis_parameters, num_sims):
        start = time.perf_counter()

        Env = Environment(
            date=(2024,7, 2, 12),
            latitude=df.loc["latitude"].iloc[1],
            longitude=df.loc["longitude"].iloc[1],
            max_expected_height=12000,
            elevation=20, #usikker på denne 
        )
        Env.set_atmospheric_model(
            type="Reanalysis",
            file="inputs/MC_env.nc",
            dictionary="ECMWF"
        )
        print("Load env: ", time.perf_counter() - start)

        N2O_T = setting["ox_temp"] # C, target with mechanical relief valve
        fuel_T = 20
        ambient_T = 20
        K_to_C = lambda C: C-273.15

        fuel_mix = pyfluids.Mixture(
            [pyfluids.FluidsList.Water, pyfluids.FluidsList.Ethanol], [25,75 ]).with_state(
                pyfluids.Input.pressure(df.loc["fuel_pressure"].iloc[1]*1e5), # 3800000 Pa tank pressure, 38 
                pyfluids.Input.temperature(df.loc["fuel_temp"].iloc[1]))

        N2O = pyfluids.Fluid(pyfluids.FluidsList.NitrousOxide).with_state(
            pyfluids.Input.pressure(df.loc["ox_pressure"].iloc[1]*1e5), #must be over 18 bar or else it will not be liquid @ -20c
            pyfluids.Input.temperature(df.loc["ox_temp"].iloc[1])) #c 

        n2i = pyfluids.Fluid(pyfluids.FluidsList.Nitrogen).with_state(
            pyfluids.Input.pressure(df.loc["n2_pressure"].iloc[1]*1e5), # 285 bar
            pyfluids.Input.temperature(ambient_T)) # Initial nitrogen state

        n2f = pyfluids.Fluid(pyfluids.FluidsList.Nitrogen).with_state(
            pyfluids.Input.pressure(df.loc["n2_pressure"].iloc[1]*1e5), 
            pyfluids.Input.temperature(fuel_T)) # Nitrogen state after entering the other tanks

        N2_tank_length = setting["n2_tank_length"]
        N2_tank_radius = setting["n2_tank_radius"]


        fuel_tank_length  = setting["tank_length"]
        fuel_tank_radius = setting["fuel_tank_radius"]

        ox_tank_length = setting["ox_tank_length"]#m
        ox_tank_radius = setting["ox_tank_radius"]#m

        fuel_tank_geometry = CylindricalTank(fuel_tank_radius, fuel_tank_length,spherical_caps=False)
        ox_tank_geometry = CylindricalTank(ox_tank_radius, ox_tank_length, spherical_caps=False)
        n2_tank_geometry = CylindricalTank(N2_tank_radius, N2_tank_length, spherical_caps=True) # what happens here?


        N2O = Fluid(name="N2O", density=N2O.density)
        fuel = Fluid(name="liq_eth90", density=fuel_mix.density)
        gas_N2i = Fluid(name="gas_N2_initial", density=n2i.density)
        gas_N2f = Fluid(name="gas_N2_final", density=n2f.density)

        ox_mass = setting["ox_mass"]
        #ox_volume = (ox_mass*0.001)/N2O.density  #L

        prop_mass = setting["fuel_mass"]
        n2_volume = setting["n2_volume"]
        n2_mass = (n2_volume*0.001)*n2i.density #kg 0.001 er for å konvertere til m^3


        of = setting["of"] # oxidizer to fuel ratio
        ox_perc = of / (1+of) # Oxidizer percentage
        prop_perc = 1 - ox_perc # Fuel percentage

        mdot=setting["mdot"]#kg/s

        ox_mdot = (mdot*ox_perc) 
        prop_mdot = mdot*prop_perc 

        if  (prop_mass/prop_mdot)*0.99999 <  (ox_mass/ox_mdot)*0.99999:
            burnout_time = (prop_mass/prop_mdot)*0.99999
        else:
            burnout_time = (ox_mass/ox_mdot)*0.99999

        print(burnout_time)
        n2_mdot = n2_mass / burnout_time

        #burnout_time /=2 ##for half burn

        ox_tank_vol = ox_tank_geometry.total_volume
        N2O_tank = MassFlowRateBasedTank(
            name="oxidizer tank",
            geometry=ox_tank_geometry,
            liquid=N2O,
            gas=gas_N2f,
            flux_time=burnout_time,

            initial_liquid_mass=ox_mass,
            initial_gas_mass=0,
            liquid_mass_flow_rate_in=0,
            liquid_mass_flow_rate_out=lambda t: ox_mdot if t < burnout_time else 0,
            gas_mass_flow_rate_in=lambda t: ox_perc*n2_mdot/2 if t < burnout_time else 0,
            gas_mass_flow_rate_out=0,
        )


        fuel_tank_volume = fuel_tank_geometry.total_volume
        fuel_tank = MassFlowRateBasedTank(
            name="fuel tank",
            geometry=fuel_tank_geometry,
            liquid=fuel,
            gas=n2f,
            flux_time=burnout_time,

            initial_liquid_mass=prop_mass,
            initial_gas_mass=0,
            liquid_mass_flow_rate_in=0,
            liquid_mass_flow_rate_out=lambda t: prop_mdot if t < burnout_time else 0,
            gas_mass_flow_rate_in=lambda t: prop_perc*n2_mdot/2 if t < burnout_time else 0,
            gas_mass_flow_rate_out=0,
        )

        n2_tank_volume = 9.8e-3 # m^3 (9.8L)
        n2_tank = MassFlowRateBasedTank(
            name="N2 tank",
            geometry=n2_tank_geometry,
            liquid=gas_N2i,
            gas=gas_N2i,
            flux_time=burnout_time,

            initial_liquid_mass=0,
            initial_gas_mass=n2_mass,
            liquid_mass_flow_rate_in=0,
            liquid_mass_flow_rate_out=0,
            gas_mass_flow_rate_in=0,
            gas_mass_flow_rate_out=lambda t: n2_mdot if t < burnout_time else 0,
        )

        print("Load tanks: ", time.perf_counter() - start)

        thrust = pd.read_csv("inputs/rocketpyeng.csv")
        thrust.iat[2,0] = burnout_time-0.5
        thrust.iat[3,0] = burnout_time


        thrust.to_csv("inputs/rocketpyeng.csv", index = False)

        print("Load thrust: ", time.perf_counter() - start)

        liquid_motor = LiquidMotor(
            thrust_source=r"inputs\rocketpyeng.csv",
            center_of_dry_mass_position=0,
            dry_inertia=(0, 0, 0),
            dry_mass=0.1,
            burn_time=(0, burnout_time),
            nozzle_radius=0.068,
            nozzle_position=0,
            coordinate_system_orientation="nozzle_to_combustion_chamber",
        )

        liquid_motor.add_tank(fuel_tank, position=setting["fuel_tank_position"]) # 
        liquid_motor.add_tank(N2O_tank, position=setting["ox_tank_position"])
        liquid_motor.add_tank(n2_tank, position=setting["n2_tank_position"])

        heimdal = Rocket(
            radius = setting["radius"],
            mass = setting["dry_mass"],
            inertia=(setting["Ixx"], setting["Iyy"], setting["Izz"], setting["Ixz"], setting["Ixy"], setting["Iyz"]),
            power_off_drag=r"inputs/drag.off.csv",
            power_on_drag=r"inputs/drag.on.csv",
            center_of_mass_without_motor=setting["cg"],
            coordinate_system_orientation="nose_to_tail"
        )

        print("Load rocket: ", time.perf_counter() - start)

        NoseCone = heimdal.add_nose(
            length = setting["nose_length"],
            kind ="von karman",
            position = 0  
        )

        heimdal.set_rail_buttons(
            upper_button_position = setting["rocket_length"]-3.1,
            lower_button_position = setting["rocket_length"]-1.1
        )


        fin_set = heimdal.add_trapezoidal_fins(
            n=4,
            root_chord=setting["fin_root_chord"],
            tip_chord=setting["fin_tip_chord"],
            span=setting["fin_span"],
            position=setting["fin_position"],
            cant_angle=0,
            sweep_length=setting["sweep_length"]
        )

        heimdal.add_motor(liquid_motor, position = setting["rocket_length"])

        drogue_chute = heimdal.add_parachute(
            'Drogue',
            cd_s=setting["drogue_Cd_s"],
            trigger=drogue_trigger,
            sampling_rate=105,
            lag=setting["drogue_lag"]
        )

        main_chute = heimdal.add_parachute(
            'Main',
            cd_s=setting["main_Cd_s"],
            trigger=main_trigger,
            sampling_rate=105,
            lag=setting["main_lag"]
        )

        print("Load rocket components: ", time.perf_counter() - start)
            

        try: 
            test_flight = Flight(
                rocket = heimdal,
                environment = Env,
                rail_length = setting["rail_length"],
                inclination = setting["inclination"],
                heading = setting["heading"],
                max_time = 1500,
                terminate_on_apogee = False,
            )

            print("Run sim: ", time.perf_counter() - start)
            print("x: ", test_flight.x_impact, "y: ", test_flight.y_impact)
            MC_sim_result.append(test_flight)
        except Exception as E:
            print("Error during simulation:", E)
            continue

    return MC_sim_result