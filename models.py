from dataclasses import dataclass

@dataclass(slots=True)
class MCParams:
    num: int

@dataclass(slots=True)
class Flight:
    x_impact: float
    y_impact: float