# Multi-Years Spawners
Fit spawner abundance and timing for multiple years, from multiple observations per year.
## Objective:
Fit a spawner abundance model to one spawning site using many years of data at once. For example, wild Sockeye spawning in the Okanagan River above Osoyoos Lake.

## Model:
Simulates spawner population dynamics including productivity (smolts/spawner), age compostion of smolts and adults, and smolt survival by ocean entry year. Within a year, spawners enter according to annually varying entry time and spread, and stay for a fixed residence time. Spawners are observed with sampling error, a random number of time each year, with random sampling dates each year. 

Fits a Stan model: a multi-level regression where level 1 parameters vary by year and share information through a distribution fitted to each level 1 parameters. That distribution is the level 2 parameters and the basis for sampling the probability density distributions (fitting) at level 1. Level 2 is the overall pattern across years and level 1 is the change from that each year. Information is shared between years, that improves the precision of parameters estimates for each year by using information from the other years; estimates for poorly observed years will tend toward the mean run for all years. 
