# uguca

**Spectral Boundary Integral Method** to simulate rupture propagation along a weak interface between two semi-infinite half-planes.

## License

Copyright &copy; 2021 ETH Zurich (David S. Kammer)

uguca is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

uguca is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with uguca.  If not, see <https://www.gnu.org/licenses/>.


## Getting started

**Requirements**:

- [CMake](https://cmake.org/) (3.1.0 or higher)
- [FFTW](http://www.fftw.org) (3.x)
- [OpenMPI](https://www.open-mpi.org/)
- [GNU Scientific Library](https://www.gnu.org/software/gsl/)


**Compilation**:

To compile uguca, follow these steps

1. `git clone https://gitlab.com/uguca/uguca.git -b stable`
2. `cd uguca`
3. `mkdir build; cd build;`
4. `cmake -DCMAKE_BUILD_TYPE:STRING=Release ..`
5. `make`


**Test**:

You may test if uguca is running well on your computer by following these steps

1. `cd build`
2. `ctest`


**Example Simulation**:

You may run an example simulation by following these steps

1. `cd build/examples`
2. `make`
3. _e.g._ `./basic_example`
4. visualize result (you need Python) _e.g._ `./basic_example_plot.py`


## Documentation

Detailed information can be found in the [online documentation](https://uguca.gitlab.io/uguca/)