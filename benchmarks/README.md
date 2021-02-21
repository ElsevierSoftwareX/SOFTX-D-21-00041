# SCEC Benchmarks

In the part, we provide example simulations that correspond to some of the benchmark simulation of the [The SCEC/USGS Spontaneous Rupture Code Verification Project](http://scecdata.usc.edu/cvws/).

To run and verify the results, go into the folder of your TPV example of choice and follow these steps:

First, **execute** the benchmark problem. Either, run directly the executable by `./TPV<id>` where `<id>` is the number of the example, or you may use the provided shell script `./TPV<id>.sh`, which uses `MPI` to run the benchmark simulation. 

Some notes:

- you may adapt various simulation parameters via command-line options. To see available options, use `./TPV<id> -h`.
- the benchmark simulations use by default a fine discretization for comparison with available references data. Hence, the simulations require considerable amount of memory. You may reduce the discretization via command-line options.

Further, you can **analyze** the simulation results using the provided `./inspect_results.py <bname>`. The `<bname>` is generated automatically by the simulation. Look for the `<bname>.info` file to determine the name of your simulation.

Finally, you may **compare** the simulation results with existing benchmark data. You may download the available data by goint into the ref folder via `cd ref` and execute `./crawler.py` script. If the data are downloaded, the  `./inspect_results.py <bname>` will automatically include the reference data.
