from dataclasses import dataclass
from datetime import datetime

@dataclass(slots=True)
class MCParams:
    num: int
    time: str

@dataclass(slots=True)
class Flight:
    x_impact: float
    y_impact: float

@dataclass(slots=True)
class RocketParams:
    elevation: tuple[float, float] = (20., 1.)  # m

@dataclass(slots=True)
class Weather:
    wind_speed: float
    wind_direction: float
    temperature: float
    pressure: float