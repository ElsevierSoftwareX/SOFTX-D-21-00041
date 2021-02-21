User's Guide
============

Example simulation
------------------

The following is the main file for the `basic_example`, which includes all central parts of uguca.
     
.. literalinclude:: ../../../examples/basic_example.cc
   :language: cpp


Key components of every simulation
----------------------------------

As you can see from the example above, a typical simulation contains the following parts:

**Material**

The material is linear elastic and defined by the elastic modulus `E`, Poisson's ratio `nu`, and density `rho`, which are provided as arguments to the constructor of the ``Material`` object. Additionally, for 2D problems one can either assume the `plane strain` (default) or `plane stress` hypothesis.

The material needs to read the pre-computed kernels required for the spectral-boundary-integral method. If the pre-computed kernels do not exist already for your choice of material properties, you can generate them (see below).

**Mesh**

The discretization of the simulation domain is given by the ``Mesh`` object. The discretization needs to be a regular mesh with constant element size due to the spectral nature of the method. The mesh is given by the length and number of elements. A 2D problem consists of a 1D interface, and, hence, is described by single length. A 3D problem, however, results in a 2D interface, which is defined by the length and width, as well as, number of elements for each side (see other examples for 3D problems).

The ``Mesh`` conveniently generated the mesh and discretization itself. However, if needed, one can also provide a pre-defined mesh.

**Interface Law**

The interface law is the constitutive law describing the interactions across the interface. It takes as arguments the mesh and the parameters for the constitutive law, *e.g.*, fracture energy and peak strength.

Various constitutive laws are available, including

- Linear Shear Cohesive Law
- Linear Coulomb Friction Law
- Mixed-mode Fracture Law (called Barras Law)
- Rate- and State- Friction Law

and correspond to the class of the law object. Detailed descriptions of the available interface laws are provided in :doc:`interface_laws`

Additional interface constitutive laws can be implemented, see :doc:`developer_guide`.

**Interface**

The interface is the core class that combines all of the above objects. Hence, its arguments include the `Mesh`, `Material`, and `Law`. Various types of interfaces exist, which include:

- ``UnimatShearInterface``, which is a shear interface between materials of same properties
- ``BimatInterface``, which is an interface between two materials of different properties
- ``DefRigInterface``, which is an interface between a deformable material and a rigid flat surface

Additional interfaces can be implmented, see :doc:`developer_guide`.

The interface needs to be initialized with the ``init()`` method after the time step is set (see below).

**External Loading**

External interface load can be applied by accessing the normal and shear load arrays of the interface using the ``getNormalLoad()`` and/or the ``getShearLoad()`` methods of the ``Interface`` object. If non-uniform profiles of loading are needed, we can use the ``getCoords()`` method of the ``Mesh`` object for coordinates of the interface discretization.

**Time Step**

The interface requires a time step used for time integration. It can be set via the ``setTimeStep(time_step)`` method of the interface. The imposed time step needs to be considerably smaller than the stable time step, which can be accessed via the ``getStableTimeStep()`` method. For 2D simulations (*i.e.*, 1D interface), it is common to use a time step that is 40% of the stable time step.

**Dumping**

The interface does dump various interface fields, *e.g.*, interface cohesion, displacement, velocity, etc. It first needs to be initialized with a `name` and `path`. The fields of interest need then to be registered and anytime a dump should be made, the ``dump(step,time)`` method of the ``Interface`` object should be called.

**Time Integration**

Finally, time integration is done with the ``advanceTimeStep()`` method of the ``Interface`` object. 


Precompute kernels for material properties
------------------------------------------

If you use material properties for which no pre-computed kernel exist in uguca, you may compute and save a new pre-computed kernel for your material properties of choice. The ``kernels/README.md`` provides additional information. You may follow these steps to generate new kernels (shown here with example parameters)::

  > cd kernels
  > ./laplace_inversion.py 0.4 pstrain 100

The first argument corresponds to the Poisson's ratio, the second determines if you use plane strain (``pstrain``) or plane stress (``pstress``) assumptions, and the last argument is the cut-off length. It is standard to use 100 for the cut-off, which is also the default value for uguca. However, you may increase it, if required. Some kernels with 200 and 300 cut-off are already provided in the ``uguca/kernels/tcut200`` and ``uguca/kernels/tcut300`` folders, respectively. If you need to use these longer kernels, you need to provide the path to these folders to the ``readPrecomputedKernels()`` method of the ``Material`` object.

We advice that you verify the kernel using the ``kernel_plotter.py`` script with arguments *Poisson ratio* and *pstrain/pstress*. Please check that there are not any random zeros in the kernel, which would indicate that the algorithm did not converge and that you need to increase the number of iterations, *i.e.*, ``mp.dps``.

Note that the ``laplace_inversion.py`` script is very slow. However, remember that you only need to do this once for your parameter choice.

	      
Develop your own simulation
---------------------------

If you wish to develop your own simulation, you need to follow these steps:

1. write your main file, *e.g.*, ``uguca/examples/my_example.cc``. 
2. add your main file to the ``uguca/examples/CMakeLists.txt``: *e.g.*, ``add_simulation(my_example my_example.cc)``
3. recompile uguca::

     > cd build
     > make

4. run your simulation::

     > cd build/examples
     > ./my_example
