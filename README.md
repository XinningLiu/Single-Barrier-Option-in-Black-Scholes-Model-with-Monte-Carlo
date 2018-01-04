# Single-Barrier-Option-in-Black-Scholes-Model-with-Monte-Carlo
  This project simulated a single barrier option (up and out option)price using Monte-Carlo method in R with variance reduction method and analyzed the difference between the theoretical price and simulated price.
  We did basic Monte-Carlo simulation under Euler Scheme, also we did two variance reduction techniqe, rescaled by value and rescaled by weight.
  Note that this basic Monte-Carlo simulation is a biased method. In case of the unbiased simulation, we need to add Brownian Bridge into the model.

Draw Basic Monte Carlo.R: This file simulated an up-and-out option using basic Monte Carlo method, and graph the result.
Comparison of 3 methods.R: Includes one basic Monte-Carlo simulation and two variance reduction simulation.
