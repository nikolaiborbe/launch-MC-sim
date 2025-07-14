import datetime
from rocketpy import Flight  # type: ignore
from fastapi import FastAPI, Depends
from pydantic import BaseModel, Field
from typing import List

import models
from weather import get_weather
from utils import filter_flight_data
from mc import get_mc_sim_result


class MCParams(BaseModel):
    num_sims: int = Field(
        10, gt=0, description="Number of Monte Carlo simulations to run."
    )
    launch_time: datetime.datetime = Field(
        default_factory=lambda: datetime.datetime.now(datetime.timezone.utc),
        description="Time of the rocket launch.",
    )
    rocket_params: models.RocketParams = Field(
        default_factory=models.RocketParams,
        description="Parameters for the rocket being simulated.",
    )


app = FastAPI(
    title="MC Simulation Server",
    description="Simulates rocket launches using Monte Carlo methods.",
)


@app.post("/mc-sim", response_model=List[models.FlightData])
async def get_mc_sim(mc_params: MCParams):

    weather = await get_weather(mc_params.launch_time)

    result: list[Flight] = get_mc_sim_result(
        mc_params.rocket_params,
        weather,
        mc_params.launch_time,
    )

    data: List[models.FlightData] = [filter_flight_data(flight) for flight in result]

    return data


@app.get("/")
async def root():
    return {
        "message": "Server for Monte Carlo simulations of rocket launches.",
    }
