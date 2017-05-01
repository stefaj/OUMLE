# OUMLE
## Maximum Likelihood Parameter Estimation for Ornstein Uhlenbeck

Based on Maximum likelihood estimation of mean reverting processes by Jose Carlos Garcia Franco

Solves for mu, theta and sigma in the equation `mu(theta - x) dt + sigma dB`.

OUMLE.Solve to solve for the parameters

```
Solve(double[] xs, double delta, out double mu, out double theta, out double sigma)
```

where `xs` is the input array, `delta` is the time different `t1 -t0`, `mu` is the output mean reversion parameters, `theta` is the output mean parameters and `sigma` is the output variance

TODO: Improve the FindEta function to use a Quasi Newton method to maximize `V(d)`

Works for now
