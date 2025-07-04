from urllib import response
from rocketpy import Flight
from fastapi import FastAPI,  Depends
from utils import filter_flight_data
import models
from MC import get_MC_sim_result

app = FastAPI(
    title="MC Simulation Server",
    description="Simulates rocket launches using Monte Carlo methods.",
)


@app.get("/MC-sim")
def get_MC_sim(MC_params: models.MCParams = Depends()):
    # list of flight sims
    result: list[Flight] = get_MC_sim_result(MC_params.num)
    # only return the "valueable" data from the flight sims

    data: list[models.Flight] = [filter_flight_data(flight) for flight in result]

    return {"data": data}


@app.get("/")
def root():
    return {
        "message": "Server for Monte Carlo simulations of rocket launches.",
    }
