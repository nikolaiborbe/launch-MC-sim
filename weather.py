import requests
from datetime import datetime
from zoneinfo import ZoneInfo


from models import Weather


USER_AGENT = "my-cool-app/1.0 " "(contact: nikolaiborbe@gmail.com)"


def _ts_to_weather(ts: dict) -> Weather:
    inst = ts["data"]["instant"]["details"]
    return Weather(
        wind_speed=inst["wind_speed"],
        wind_direction=inst["wind_from_direction"],
        temperature=inst["air_temperature"],
        pressure=inst["air_pressure_at_sea_level"],
    )


async def get_weather(t: datetime, lat: float = 63.786667, lon: float = 9.363056) -> Weather:
    """
    Return the first forecast *at or after* the UTC instant *t*
    for the coordinates defined by lat, lon.

    `t` may be naive (assumed local Europe/Oslo) or timezone-aware.
    Raises `ValueError` if the requested time is outside the served range.
    """
    if t.tzinfo is None:
        t = t.replace(tzinfo=ZoneInfo("Europe/Oslo"))
    t_utc = t.astimezone(ZoneInfo("UTC"))

    url = (
        "https://api.met.no/weatherapi/locationforecast/2.0/compact"
        f"?lat={lat}&lon={lon}"
    )
    r = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=10)
    r.raise_for_status()
    timeseries = r.json()["properties"]["timeseries"]

    def ts_dt(item):  # → aware UTC datetime
        return datetime.fromisoformat(item["time"].replace("Z", "+00:00"))

    try:
        ts = next(ts for ts in sorted(timeseries, key=ts_dt) if ts_dt(ts) >= t_utc)
    except StopIteration:
        raise ValueError("Requested time is outside the available forecast window")

    return _ts_to_weather(ts)
