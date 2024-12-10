# Multilevel-Spawners
## Objective:
Analyze surveys of spawning salmon: fit spawner abundance, timing, spread, and residence.
## Case
Wild Sockeye Salmon spawning in the Okanagan River above Osoyoos Lake (natal lake); 9 to 19 surveys/year, 2000 to 2023, total 227.
## Models:
**Simulations:** spawner population dynamics including productivity (smolts/spawner), age compostion of smolts and spawners, and smolt survival by ocean entry year. Within a year, spawners enter according to annually varying entry timing and spread, and stay for a fixed residence time. Spawners are observed with sampling error, a random number of time each year, with random sampling dates each year.
**Stan models:** non-linear multi-level regression (Stan) where parameters vary by year but share information through a fitted distribution of parameters across years. That distribution is the basis for sampling annual PDDs of parameters (fitting year values). Sharing  information between years in this way improves precision and location of parameter estimates for poorly observed years.
