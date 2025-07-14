import requests

payload = {
    "num_sims": 10,
    "rocket_params": {
        "mass": 12_000,
        "burn_time": 180,
        "thrust": 2.5e6,
    },
}

resp = requests.post("http://localhost:8000/mc-sim", json=payload, timeout=15)
resp.raise_for_status()
print(resp.json())
