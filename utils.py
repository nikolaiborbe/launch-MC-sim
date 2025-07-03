from rocketpy import Flight
from typing import Union

def drogue_trigger(p, h, y):
    # Check if rocket is going down, i.e. if it has passed the apogee
    vertical_velocity = y[5]
    # Return true to activate parachute once the vertical velocity is negative
    return vertical_velocity < 0 # False if you want ballistic

def main_trigger(p,h,y):
    vertical_velocity = y[5]
    return y[5] < 0 and h <= 400 # False if you want ballistic

def filter_flight_data(flight: Flight) -> dict:
    return {
        "velocity": flight.apogee,
    }