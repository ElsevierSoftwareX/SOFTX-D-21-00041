Getting Started
===============

Obtaining uguca
---------------

The easiest way of obtaining uguca is to clone it from `gitlab <https://gitlab.com/>`_::

  > git clone git@gitlab.com:uguca/uguca.git -b stable

If you do not have a gitlab account, you should use::

  > git clone https://gitlab.com/uguca/uguca.git -b stable

Alternatively, if you are not familiar with ``git``, you may download the source code of stable releases at: `uguca release <https://gitlab.com/uguca/uguca/-/releases>`_


Requirements for uguca
----------------------

The following software are required for uguca:

- `CMake <https://cmake.org/>`_ (3.1.0 or higher)
- `FFTW <http://www.fftw.org>`_ (3.x)
- `OpenMPI <https://www.open-mpi.org/>`_
- `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_
  
Optional software for additional features in uguca:

- `git <https://git-scm.com/>`_
- `Automatically Tuned Linear Algebra Software (ATLAS) <http://math-atlas.sourceforge.net/>`_
- `Python3 <https://www.python.org/>`_
- `CMake curses graphical user interface <https://cmake.org/>`_

Compilers:

uguca has successfully been compiled with the following compilers:

- gcc-10.2
- gcc-7.5
- gcc-5.4
- clang-12
  

**Ubuntu-based systems**:

In ubuntu-based systems you can install the requirements with the following command::

  > sudo apt-get install g++ cmake git libfftw3-dev libopenmpi-dev openmpi-bin libgsl-dev 

Uguca includes additional features, which have further requirements.

Required packages for ccmake::

  > sudo apt-get install cmake-curses-gui

Required packages for UCA_USE_BLAS option::

  > sudo apt-get install libatlas-base-dev

Required packages for doc option (doc is available online)::

  > sudo apt-get install python3 python3-sphinx python3-sphinx-rtd-theme

If standard FFTW installation does not work refer to   :doc:`this guide <FFTW>`.
  
**macOS**:

First, install Command Line Tools::

  > xcode-select --install

You can install the requirements with `Homebrew <https://brew.sh>`_::

  > brew install cmake fftw gsl 

Required packages for *UCA_USE_MPI* option::

  > brew install open-mpi

Required packages for doc option (doc is online available)::

  > brew install python3
  > pip3 install sphinx sphinx-rtd-theme

Compiling uguca
---------------

**Ubuntu-based and macOS systems**:

You can configure and build uguca by following these steps::

  > cd uguca
  > mkdir build
  > cd build
  > cmake -DCMAKE_BUILD_TYPE:STRING=Release ..
  > make

If you would like to run the uguca tests as verification::

  > ctest

  
Running an example
------------------

Several example simulations are provided in the `uguca/examples` folder. To run a simulation, you typically proceed as follows::

  > cd build/examples
  > make
  > ./basic_example
  
You may visualize the results with the provided script::

  > ./basic_example_plot.py

The process of developing and running your own simulations is described in details in :doc:`./user_guide`.
