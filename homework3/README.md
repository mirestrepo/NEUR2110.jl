## Homework 3: Multivariate processes; Power, Coherence and Partial Coherence

This homework consists of two main problems

The plots require javascript rendering, which GitHub does not allow. Therefore for proper rendering of the solutions, please referred to the [nbviewer version](http://nbviewer.jupyter.org/github/mirestrepo/NEUR2110.jl/blob/master/homework1):

* [Problem 1 and 2](http://nbviewer.jupyter.org/github/mirestrepo/NEUR2110.jl/blob/master/homework3/Problem1and2.ipynb): VAR(3), Stability, Coherence and Partial Coherence Power Based on Multitaper Spectra and Maximum Entropy Estimated VAR(p) Model

Note: Three additional files are included

1. multitaper.jl: Perform tapered FFT of continuous signals - Originally writteb by [Simon Kornblith](https://github.com/simonster/). It is an updated version from this [here](https://github.com/simonster/Synchrony.jl/blob/403bf6d4e65bab45b5806787e8efa48e390a86f7/src/multitaper.jl)
2. transform_stats.jl: Tools for statistical computations between sets of signals. Originally writteb by [Simon Kornblith](https://github.com/simonster/). It is an updated version from this [here](https://github.com/simonster/Synchrony.jl/blob/403bf6d4e65bab45b5806787e8efa48e390a86f7/src/transform_stats.jl)
3. var\_maxent.jl: Fits Vector AR model (VAR) to multi-trial, multivariate time series data using maximum entropy estimation. This file is basically a direct translation from the file provided in class var\_maxent.m. Note that the code does not take advantage of Julia specific declarations and may be inefficient.


