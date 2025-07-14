from dataclasses import dataclass


@dataclass(slots=True)
class FlightData:
    x_impact: float
    y_impact: float


@dataclass(slots=True)
class RocketParams:
    elevation: float = 20.0  # m


@dataclass(slots=True)
class Weather:
    wind_speed: float
    wind_direction: float
    temperature: float
    pressure: float
