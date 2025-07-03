from urllib import response
from rocketpy import Flight
from fastapi import FastAPI, HTTPException, Depends 
from utils import filter_flight_data
from contextlib import asynccontextmanager
from models import MCParams
from MC import get_MC_sim_result

app = FastAPI(
    title="MC Simulation Server",
    description="Simulates rocket launches using Monte Carlo methods.",
)

@app.get("/MC-sim")
def get_MC_sim(MC_params: MCParams = Depends()):
    # list of flight sims
    result: list[Flight] = get_MC_sim_result(MC_params.num)
    # only return the "valueable" data from the flight sims

    data = [filter_flight_data(flight) for flight in result]

    return {
        "message": "Rocket simulation API. Use GET /MC-sim to fetch current metrics.",
        "params": MC_params,
        "data": data
    }


