Package details
===============

Here is a brief description of all the standard and user packages in
LAMMPS.  It lists authors (if applicable) and summarizes the package
contents.  It has specific instructions on how to install the package,
including, if necessary, info on how to download or build any extra
library it requires.  It also gives links to documentation, example
scripts, and pictures/movies (if available) that illustrate use of the
package.

The majority of packages can be included in a LAMMPS build with a
single setting (``-D PGK_<NAME>=on`` for CMake) or command
(``make yes-<name>`` for make).  See the :doc:`Build package <Build_package>`
page for more info.  A few packages may require additional steps;
this is indicated in the descriptions below.  The :doc:`Build extras <Build_extras>`
page gives those details.

.. note::

   To see the complete list of commands a package adds to LAMMPS,
   you can examine the files in its src directory, e.g. "ls
   src/GRANULAR".  Files with names that start with fix, compute, atom,
   pair, bond, angle, etc correspond to commands with the same style name
   as contained in the file name.

.. table_from_list::
   :columns: 6

   * :ref:`ASPHERE <PKG-ASPHERE>`
   * :ref:`BODY <PKG-BODY>`
   * :ref:`CLASS2 <PKG-CLASS2>`
   * :ref:`COLLOID <PKG-COLLOID>`
   * :ref:`COMPRESS <PKG-COMPRESS>`
   * :ref:`CORESHELL <PKG-CORESHELL>`
   * :ref:`DIPOLE <PKG-DIPOLE>`
   * :ref:`GPU <PKG-GPU>`
   * :ref:`GRANULAR <PKG-GRANULAR>`
   * :ref:`KIM <PKG-KIM>`
   * :ref:`KOKKOS <PKG-KOKKOS>`
   * :ref:`KSPACE <PKG-KSPACE>`
   * :ref:`LATTE <PKG-LATTE>`
   * :ref:`MANYBODY <PKG-MANYBODY>`
   * :ref:`MC <PKG-MC>`
   * :ref:`MESSAGE <PKG-MESSAGE>`
   * :ref:`MISC <PKG-MISC>`
   * :ref:`MLIAP <PKG-MLIAP>`
   * :ref:`MOLECULE <PKG-MOLECULE>`
   * :ref:`MPIIO <PKG-MPIIO>`
   * :ref:`MSCG <PKG-MSCG>`
   * :ref:`OPT <PKG-OPT>`
   * :ref:`PERI <PKG-PERI>`
   * :ref:`PLUGIN <PKG-PLUGIN>`
   * :ref:`POEMS <PKG-POEMS>`
   * :ref:`PYTHON <PKG-PYTHON>`
   * :ref:`QEQ <PKG-QEQ>`
   * :ref:`REPLICA <PKG-REPLICA>`
   * :ref:`RIGID <PKG-RIGID>`
   * :ref:`SHOCK <PKG-SHOCK>`
   * :ref:`SNAP <PKG-SNAP>`
   * :ref:`SPIN <PKG-SPIN>`
   * :ref:`SRD <PKG-SRD>`
   * :ref:`VORONOI <PKG-VORONOI>`

.. table_from_list::
   :columns: 6

   * :ref:`USER-ADIOS <PKG-USER-ADIOS>`
   * :ref:`USER-ATC <PKG-USER-ATC>`
   * :ref:`USER-AWPMD <PKG-USER-AWPMD>`
   * :ref:`USER-BOCS <PKG-USER-BOCS>`
   * :ref:`USER-CGDNA <PKG-USER-CGDNA>`
   * :ref:`USER-CGSDK <PKG-USER-CGSDK>`
   * :ref:`USER-COLVARS <PKG-USER-COLVARS>`
   * :ref:`USER-DIFFRACTION <PKG-USER-DIFFRACTION>`
   * :ref:`USER-DPD <PKG-USER-DPD>`
   * :ref:`USER-DRUDE <PKG-USER-DRUDE>`
   * :ref:`USER-EFF <PKG-USER-EFF>`
   * :ref:`USER-FEP <PKG-USER-FEP>`
   * :ref:`USER-H5MD <PKG-USER-H5MD>`
   * :ref:`USER-INTEL <PKG-USER-INTEL>`
   * :ref:`USER-LB <PKG-USER-LB>`
   * :ref:`USER-MANIFOLD <PKG-USER-MANIFOLD>`
   * :ref:`USER-MEAMC <PKG-USER-MEAMC>`
   * :ref:`USER-MESODPD <PKG-USER-MESODPD>`
   * :ref:`USER-MESONT <PKG-USER-MESONT>`
   * :ref:`USER-MGPT <PKG-USER-MGPT>`
   * :ref:`USER-MISC <PKG-USER-MISC>`
   * :ref:`USER-MOFFF <PKG-USER-MOFFF>`
   * :ref:`USER-MOLFILE <PKG-USER-MOLFILE>`
   * :ref:`USER-NETCDF <PKG-USER-NETCDF>`
   * :ref:`USER-OMP <PKG-USER-OMP>`
   * :ref:`USER-PACE <PKG-USER-PACE>`
   * :ref:`USER-PHONON <PKG-USER-PHONON>`
   * :ref:`USER-PLUMED <PKG-USER-PLUMED>`
   * :ref:`USER-PTM <PKG-USER-PTM>`
   * :ref:`USER-QMMM <PKG-USER-QMMM>`
   * :ref:`USER-QTB <PKG-USER-QTB>`
   * :ref:`USER-QUIP <PKG-USER-QUIP>`
   * :ref:`USER-REACTION <PKG-USER-REACTION>`
   * :ref:`USER-REAXC <PKG-USER-REAXC>`
   * :ref:`USER-SCAFACOS <PKG-USER-SCAFACOS>`
   * :ref:`USER-SDPD <PKG-USER-SDPD>`
   * :ref:`USER-SMD <PKG-USER-SMD>`
   * :ref:`USER-SMTBQ <PKG-USER-SMTBQ>`
   * :ref:`USER-SPH <PKG-USER-SPH>`
   * :ref:`USER-TALLY <PKG-USER-TALLY>`
   * :ref:`USER-UEF <PKG-USER-UEF>`
   * :ref:`USER-VTK <PKG-USER-VTK>`
   * :ref:`USER-YAFF <PKG-USER-YAFF>`

----------

.. _PKG-ASPHERE:

ASPHERE package
---------------

**Contents:**

Computes, time-integration fixes, and pair styles for aspherical
particle models including ellipsoids, 2d lines, and 3d triangles.

**Supporting info:**

* src/ASPHERE: filenames -> commands
* :doc:`Howto spherical <Howto_spherical>`
* :doc:`pair_style gayberne <pair_gayberne>`
* :doc:`pair_style resquared <pair_resquared>`
* `doc/PDF/pair_gayberne_extra.pdf <PDF/pair_gayberne_extra.pdf>`_
* `doc/PDF/pair_resquared_extra.pdf <PDF/pair_resquared_extra.pdf>`_
* examples/ASPHERE
* examples/ellipse
* https://lammps.sandia.gov/movies.html#line
* https://lammps.sandia.gov/movies.html#tri

----------

.. _PKG-BODY:

BODY package
------------

**Contents:**

Body-style particles with internal structure.  Computes,
time-integration fixes, pair styles, as well as the body styles
themselves.  See the :doc:`Howto body <Howto_body>` page for an
overview.

**Supporting info:**

* src/BODY filenames -> commands
* :doc:`Howto_body <Howto_body>`
* :doc:`atom_style body <atom_style>`
* :doc:`fix nve/body <fix_nve_body>`
* :doc:`pair_style body/nparticle <pair_body_nparticle>`
* examples/body

----------

.. _PKG-CLASS2:

CLASS2 package
--------------

**Contents:**

Bond, angle, dihedral, improper, and pair styles for the COMPASS
CLASS2 molecular force field.

**Supporting info:**

* src/CLASS2: filenames -> commands
* :doc:`bond_style class2 <bond_class2>`
* :doc:`angle_style class2 <angle_class2>`
* :doc:`dihedral_style class2 <dihedral_class2>`
* :doc:`improper_style class2 <improper_class2>`
* :doc:`pair_style lj/class2 <pair_class2>`

----------

.. _PKG-COLLOID:

COLLOID package
---------------

**Contents:**

Coarse-grained finite-size colloidal particles.  Pair styles and fix
wall styles for colloidal interactions.  Includes the Fast Lubrication
Dynamics (FLD) method for hydrodynamic interactions, which is a
simplified approximation to Stokesian dynamics.

**Authors:** This package includes Fast Lubrication Dynamics pair styles
which were created by Amit Kumar and Michael Bybee from Jonathan
Higdon's group at UIUC.

**Supporting info:**

* src/COLLOID: filenames -> commands
* :doc:`fix wall/colloid <fix_wall>`
* :doc:`pair_style colloid <pair_colloid>`
* :doc:`pair_style yukawa/colloid <pair_yukawa_colloid>`
* :doc:`pair_style brownian <pair_brownian>`
* :doc:`pair_style lubricate <pair_lubricate>`
* :doc:`pair_style lubricateU <pair_lubricateU>`
* examples/colloid
* examples/srd

----------

.. _PKG-COMPRESS:

COMPRESS package
----------------

**Contents:**

Compressed output of dump files via the zlib compression library,
using dump styles with a "gz" in their style name.

To use this package you must have the zlib compression library
available on your system.

**Author:** Axel Kohlmeyer (Temple U).

**Install:**

This package has :ref:`specific installation instructions <compress>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/COMPRESS: filenames -> commands
* src/COMPRESS/README
* lib/compress/README
* :doc:`dump atom/gz <dump>`
* :doc:`dump cfg/gz <dump>`
* :doc:`dump custom/gz <dump>`
* :doc:`dump xyz/gz <dump>`

----------

.. _PKG-CORESHELL:

CORESHELL package
-----------------

**Contents:**

Compute and pair styles that implement the adiabatic core/shell model
for polarizability.  The pair styles augment Born, Buckingham, and
Lennard-Jones styles with core/shell capabilities.  The :doc:`compute temp/cs <compute_temp_cs>` command calculates the temperature of a
system with core/shell particles.  See the :doc:`Howto coreshell <Howto_coreshell>` page for an overview of how to use
this package.

**Author:** Hendrik Heenen (Technical U of Munich).

**Supporting info:**

* src/CORESHELL: filenames -> commands
* :doc:`Howto coreshell <Howto_coreshell>`
* :doc:`Howto polarizable <Howto_polarizable>`
* :doc:`compute temp/cs <compute_temp_cs>`
* :doc:`pair_style born/coul/long/cs <pair_cs>`
* :doc:`pair_style buck/coul/long/cs <pair_cs>`
* :doc:`pair_style lj/cut/coul/long/cs <pair_lj>`
* examples/coreshell

----------

.. _PKG-DIPOLE:

DIPOLE package
--------------

**Contents:**

An atom style and several pair styles for point dipole models with
short-range or long-range interactions.

**Supporting info:**

* src/DIPOLE: filenames -> commands
* :doc:`atom_style dipole <atom_style>`
* :doc:`pair_style lj/cut/dipole/cut <pair_dipole>`
* :doc:`pair_style lj/cut/dipole/long <pair_dipole>`
* :doc:`pair_style lj/long/dipole/long <pair_dipole>`
* examples/dipole

----------

.. _PKG-GPU:

GPU package
-----------

**Contents:**

Dozens of pair styles and a version of the PPPM long-range Coulombic
solver optimized for GPUs.  All such styles have a "gpu" as a suffix
in their style name. The GPU code can be compiled with either CUDA or
OpenCL, however the OpenCL variants are no longer actively maintained
and only the CUDA versions are regularly tested.  The
:doc:`Speed_gpu` page gives details of what hardware and GPU
software is required on your system, and details on how to build and
use this package.  Its styles can be invoked at run time via the "-sf
gpu" or "-suffix gpu" :doc:`command-line switches <Run_options>`.  See
also the :ref:`KOKKOS <PKG-KOKKOS>` package, which has GPU-enabled styles.

**Authors:** Mike Brown (Intel) while at Sandia and ORNL and Trung Nguyen
(Northwestern U) while at ORNL and later. AMD HIP support by Evgeny
Kuznetsov, Vladimir Stegailov, and Vsevolod Nikolskiy (HSE University).

**Install:**

This package has :ref:`specific installation instructions <gpu>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/GPU: filenames -> commands
* src/GPU/README
* lib/gpu/README
* :doc:`Accelerator packages <Speed_packages>`
* :doc:`GPU package <Speed_gpu>`
* :doc:`Section 2.6 -sf gpu <Run_options>`
* :doc:`Section 2.6 -pk gpu <Run_options>`
* :doc:`package gpu <package>`
* :doc:`Commands <Commands_all>` pages (:doc:`pair <Commands_pair>`, :doc:`kspace <Commands_kspace>`)
  for styles followed by (g)
* `Benchmarks page <https://lammps.sandia.gov/bench.html>`_ of web site

----------

.. _PKG-GRANULAR:

GRANULAR package
----------------

**Contents:**

Pair styles and fixes for finite-size granular particles, which
interact with each other and boundaries via frictional and dissipative
potentials.

**Supporting info:**

* src/GRANULAR: filenames -> commands
* :doc:`Howto granular <Howto_granular>`
* :doc:`fix pour <fix_pour>`
* :doc:`fix wall/gran <fix_wall_gran>`
* :doc:`pair_style gran/hooke <pair_gran>`
* :doc:`pair_style gran/hertz/history <pair_gran>`
* examples/granregion
* examples/pour
* bench/in.chute
* https://lammps.sandia.gov/pictures.html#jamming
* https://lammps.sandia.gov/movies.html#hopper
* https://lammps.sandia.gov/movies.html#dem
* https://lammps.sandia.gov/movies.html#brazil
* https://lammps.sandia.gov/movies.html#granregion

----------

.. _PKG-KIM:

KIM package
-----------

**Contents:**

This package contains a command with a set of sub-commands that serve as a
wrapper on the
`Open Knowledgebase of Interatomic Models (OpenKIM) <https://openkim.org>`_
repository of interatomic models (IMs) enabling compatible ones to be used in
LAMMPS simulations.


This includes :doc:`kim init <kim_commands>`, and
:doc:`kim interactions <kim_commands>` commands to select, initialize and
instantiate the IM, a :doc:`kim query <kim_commands>` command to perform web
queries for material property predictions of OpenKIM IMs, a
:doc:`kim param <kim_commands>` command to access KIM Model Parameters from
LAMMPS, and a :doc:`kim property <kim_commands>` command to write material
properties computed in LAMMPS to standard KIM property instance format.

Support for KIM IMs that conform to the
`KIM Application Programming Interface (API) <https://openkim.org/kim-api/>`_
is provided by the :doc:`pair_style kim <pair_kim>` command.

.. note::

   The command *pair_style kim* is called by *kim interactions* and is not
   recommended to be directly used in input scripts.

To use this package you must have the KIM API library available on your
system. The KIM API is available for download on the
`OpenKIM website <https://openkim.org/kim-api/>`_.
When installing LAMMPS from binary, the kim-api package
is a dependency that is automatically downloaded and installed.

Information about the KIM project can be found at its website:
`https://openkim.org <https://openkim.org>`_.
The KIM project is led by Ellad Tadmor and Ryan Elliott (U Minnesota)
and is funded by the `National Science Foundation <https://www.nsf.gov/>`_.

**Authors:** Ryan Elliott (U Minnesota) is the main developer for the KIM
API and the *pair_style kim* command. Yaser Afshar (U Minnesota),
Axel Kohlmeyer (Temple U), Ellad Tadmor (U Minnesota), and
Daniel Karls (U Minnesota) contributed to the
:doc:`kim command <kim_commands>` interface in close collaboration with
Ryan Elliott.

**Install:**

This package has :ref:`specific installation instructions <kim>` on the
:doc:`Build extras <Build_extras>` page.

**Supporting info:**

* :doc:`kim command <kim_commands>`
* :doc:`pair_style kim <pair_kim>`
* src/KIM: filenames -> commands
* src/KIM/README
* lib/kim/README
* examples/kim

----------

.. _PKG-KOKKOS:

KOKKOS package
--------------

**Contents:**

Dozens of atom, pair, bond, angle, dihedral, improper, fix, compute
styles adapted to compile using the Kokkos library which can convert
them to OpenMP or CUDA code so that they run efficiently on multicore
CPUs, KNLs, or GPUs.  All the styles have a "kk" as a suffix in their
style name.  The :doc:`KOKKOS package <Speed_kokkos>` page gives
details of what hardware and software is required on your system, and
how to build and use this package.  Its styles can be invoked at run
time via the "-sf kk" or "-suffix kk" :doc:`command-line switches <Run_options>`.  Also see the :ref:`GPU <PKG-GPU>`, :ref:`OPT <PKG-OPT>`,
:ref:`USER-INTEL <PKG-USER-INTEL>`, and :ref:`USER-OMP <PKG-USER-OMP>` packages, which
have styles optimized for CPUs, KNLs, and GPUs.

You must have a C++11 compatible compiler to use this package.
KOKKOS makes extensive use of advanced C++ features, which can
expose compiler bugs, especially when compiling for maximum
performance at high optimization levels. Please see the file
lib/kokkos/README for a list of compilers and their respective
platforms, that are known to work.

**Authors:** The KOKKOS package was created primarily by Christian Trott
and Stan Moore (Sandia), with contributions from other folks as well.
It uses the open-source `Kokkos library <https://github.com/kokkos>`_
which was developed by Carter Edwards, Christian Trott, and others at
Sandia, and which is included in the LAMMPS distribution in
lib/kokkos.

**Install:**

This package has :ref:`specific installation instructions <kokkos>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/KOKKOS: filenames -> commands
* src/KOKKOS/README
* lib/kokkos/README
* :doc:`Accelerator packages <Speed_packages>`
* :doc:`KOKKOS package <Speed_kokkos>`
* :doc:`Section 2.6 -k on ... <Run_options>`
* :doc:`Section 2.6 -sf kk <Run_options>`
* :doc:`Section 2.6 -pk kokkos <Run_options>`
* :doc:`package kokkos <package>`
* Search the :doc:`commands <Commands_all>` pages (:doc:`fix <Commands_fix>`, :doc:`compute <Commands_compute>`,
  :doc:`pair <Commands_pair>`, :doc:`bond, angle, dihedral, improper <Commands_bond>`,
  :doc:`kspace <Commands_kspace>`) for styles followed by (k)
* `Benchmarks page <https://lammps.sandia.gov/bench.html>`_ of web site

----------

.. _PKG-KSPACE:

KSPACE package
--------------

**Contents:**

A variety of long-range Coulombic solvers, as well as pair styles
which compute the corresponding short-range pairwise Coulombic
interactions.  These include Ewald, particle-particle particle-mesh
(PPPM), and multilevel summation method (MSM) solvers.

**Install:**

Building with this package requires a 1d FFT library be present on
your system for use by the PPPM solvers.  This can be the KISS FFT
library provided with LAMMPS, third party libraries like FFTW, or a
vendor-supplied FFT library.  See the :doc:`Build settings <Build_settings>` page for details on how to select
different FFT options for your LAMPMS build.

**Supporting info:**

* src/KSPACE: filenames -> commands
* :doc:`kspace_style <kspace_style>`
* `doc/PDF/kspace.pdf <PDF/kspace.pdf>`_
* :doc:`Howto tip3p <Howto_tip3p>`
* :doc:`Howto tip4p <Howto_tip4p>`
* :doc:`Howto spc <Howto_spc>`
* :doc:`pair_style coul <pair_coul>`
* Search the :doc:`pair style <Commands_pair>` page for styles with "long" or "msm" in name
* examples/peptide
* bench/in.rhodo

----------

.. _PKG-LATTE:

LATTE package
-------------

**Contents:**

A fix command which wraps the LATTE DFTB code, so that molecular
dynamics can be run with LAMMPS using density-functional tight-binding
quantum forces calculated by LATTE.

More information on LATTE can be found at this web site:
`https://github.com/lanl/LATTE <latte-home_>`_.  A brief technical
description is given with the :doc:`fix latte <fix_latte>` command.

.. _latte-home: https://github.com/lanl/LATTE

**Authors:** Christian Negre (LANL) and Steve Plimpton (Sandia).  LATTE
itself is developed at Los Alamos National Laboratory by Marc
Cawkwell, Anders Niklasson, and Christian Negre.

**Install:**

This package has :ref:`specific installation instructions <latte>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/LATTE: filenames -> commands
* src/LATTE/README
* lib/latte/README
* :doc:`fix latte <fix_latte>`
* examples/latte
* `LAMMPS-LATTE tutorial <https://github.com/lanl/LATTE/wiki/Using-LATTE-through-LAMMPS>`_

----------

.. _PKG-MANYBODY:

MANYBODY package
----------------

**Contents:**

A variety of many-body and bond-order potentials.  These include
(AI)REBO, BOP, EAM, EIM, Stillinger-Weber, and Tersoff potentials.

**Supporting info:**

* src/MANYBODY: filenames -> commands
* :doc:`Pair style <Commands_pair>` page
* examples/comb
* examples/eim
* examples/nb3d
* examples/shear
* examples/streitz
* examples/vashishta
* bench/in.eam

----------

.. _PKG-MC:

MC package
----------

**Contents:**

Several fixes and a pair style that have Monte Carlo (MC) or MC-like
attributes.  These include fixes for creating, breaking, and swapping
bonds, for performing atomic swaps, and performing grand-canonical MC
(GCMC) or similar processes in conjunction with dynamics.

**Supporting info:**

* src/MC: filenames -> commands
* :doc:`fix atom/swap <fix_atom_swap>`
* :doc:`fix bond/break <fix_bond_break>`
* :doc:`fix bond/create <fix_bond_create>`
* :doc:`fix bond/create/angle <fix_bond_create>`
* :doc:`fix bond/swap <fix_bond_swap>`
* :doc:`fix charge/regulation <fix_charge_regulation>`
* :doc:`fix gcmc <fix_gcmc>`
* :doc:`fix tfmc <fix_tfmc>`
* :doc:`fix widom <fix_widom>`
* :doc:`pair_style dsmc <pair_dsmc>`
* https://lammps.sandia.gov/movies.html#gcmc

----------

.. _PKG-MESSAGE:

MESSAGE package
---------------

**Contents:**

Commands to use LAMMPS as either a client or server and couple it to
another application.

**Install:**

This package has :ref:`specific installation instructions <message>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/MESSAGE: filenames -> commands
* lib/message/README
* :doc:`message <message>`
* :doc:`fix client/md <fix_client_md>`
* :doc:`server md <server_md>`
* :doc:`server mc <server_mc>`
* examples/message

----------

.. _PKG-MISC:

MISC package
------------

**Contents:**

A variety of compute, fix, pair, dump styles with specialized
capabilities that don't align with other packages.  Do a directory
listing, "ls src/MISC", to see the list of commands.

.. note::

   the MISC package contains styles that require using the
   -restrict flag, when compiling with Intel compilers.

**Supporting info:**

* src/MISC: filenames -> commands
* :doc:`compute ti <compute_ti>`
* :doc:`fix evaporate <fix_evaporate>`
* :doc:`fix orient/fcc <fix_orient>`
* :doc:`fix ttm <fix_ttm>`
* :doc:`fix thermal/conductivity <fix_thermal_conductivity>`
* :doc:`fix viscosity <fix_viscosity>`
* examples/KAPPA
* examples/VISCOSITY
* https://lammps.sandia.gov/pictures.html#ttm
* https://lammps.sandia.gov/movies.html#evaporation

----------

.. _PKG-MLIAP:

MLIAP package
-------------

**Contents:**

A general interface for machine-learning interatomic potentials, including PyTorch.

**Install:**

To use this package, also the :ref:`SNAP package <PKG-SNAP>` package needs
to be installed.  To make the *mliappy* model available, also the
:ref:`PYTHON package <PKG-PYTHON>` package needs to be installed, the version
of Python must be 3.6 or later, and the `cython <https://cython.org/>`_ software
must be installed.

**Author:** Aidan Thompson (Sandia), Nicholas Lubbers (LANL).

**Supporting info:**

* src/MLIAP: filenames -> commands
* src/MLIAP/README
* :doc:`pair_style mliap <pair_mliap>`
* :doc:`compute_style mliap <compute_mliap>`
* examples/mliap (see README)

When built with the *mliappy* model this package includes an extension for
coupling with Python models, including PyTorch. In this case, the Python
interpreter linked to LAMMPS will need the ``cython`` and ``numpy`` modules
installed.  The provided examples build models with PyTorch, which would
therefore also needs to be installed to run those examples.

----------

.. _PKG-MOLECULE:

MOLECULE package
----------------

**Contents:**

A large number of atom, pair, bond, angle, dihedral, improper styles
that are used to model molecular systems with fixed covalent bonds.
The pair styles include the Dreiding (hydrogen-bonding) and CHARMM
force fields, and a TIP4P water model.

**Supporting info:**

* src/MOLECULE: filenames -> commands
* :doc:`atom_style <atom_style>`
* :doc:`bond_style <bond_style>`
* :doc:`angle_style <angle_style>`
* :doc:`dihedral_style <dihedral_style>`
* :doc:`improper_style <improper_style>`
* :doc:`pair_style hbond/dreiding/lj <pair_hbond_dreiding>`
* :doc:`pair_style lj/charmm/coul/charmm <pair_charmm>`
* :doc:`Howto bioFF <Howto_bioFF>`
* examples/cmap
* examples/dreiding
* examples/micelle,
* examples/peptide
* bench/in.chain
* bench/in.rhodo

----------

.. _PKG-MPIIO:

MPIIO package
-------------

**Contents:**

Support for parallel output/input of dump and restart files via the
MPIIO library.  It adds :doc:`dump styles <dump>` with a "mpiio" in
their style name.  Restart files with an ".mpiio" suffix are also
written and read in parallel.

**Supporting info:**

* src/MPIIO: filenames -> commands
* :doc:`dump <dump>`
* :doc:`restart <restart>`
* :doc:`write_restart <write_restart>`
* :doc:`read_restart <read_restart>`

----------

.. _PKG-mscg:

MSCG package
------------

**Contents:**

A :doc:`fix mscg <fix_mscg>` command which can parameterize a
Multi-Scale Coarse-Graining (MSCG) model using the open-source `MS-CG library <mscg-home_>`_.

.. _mscg-home: https://github.com/uchicago-voth/MSCG-release

To use this package you must have the MS-CG library available on your
system.

**Authors:** The fix was written by Lauren Abbott (Sandia).  The MS-CG
library was developed by Jacob Wagner in Greg Voth's group at the
University of Chicago.

**Install:**

This package has :ref:`specific installation instructions <mscg>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/MSCG: filenames -> commands
* src/MSCG/README
* lib/mscg/README
* examples/mscg

----------

.. _PKG-OPT:

OPT package
-----------

**Contents:**

A handful of pair styles which are optimized for improved CPU
performance on single or multiple cores.  These include EAM, LJ,
CHARMM, and Morse potentials.  The styles have an "opt" suffix in
their style name.  The :doc:`OPT package <Speed_opt>` page gives
details of how to build and use this package.  Its styles can be
invoked at run time via the "-sf opt" or "-suffix opt" :doc:`command-line switches <Run_options>`.  See also the :ref:`KOKKOS <PKG-KOKKOS>`,
:ref:`USER-INTEL <PKG-USER-INTEL>`, and :ref:`USER-OMP <PKG-USER-OMP>` packages, which
have styles optimized for CPU performance.

**Authors:** James Fischer (High Performance Technologies), David Richie,
and Vincent Natoli (Stone Ridge Technology).

**Install:**

This package has :ref:`specific installation instructions <opt>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/OPT: filenames -> commands
* :doc:`Accelerator packages <Speed_packages>`
* :doc:`OPT package <Speed_opt>`
* :doc:`Section 2.6 -sf opt <Run_options>`
* Search the :doc:`pair style <Commands_pair>` page for styles followed by (t)
* `Benchmarks page <https://lammps.sandia.gov/bench.html>`_ of web site

----------

.. _PKG-PERI:

PERI package
------------

**Contents:**

An atom style, several pair styles which implement different
Peridynamics materials models, and several computes which calculate
diagnostics.  Peridynamics is a particle-based meshless continuum
model.

**Authors:** The original package was created by Mike Parks (Sandia).
Additional Peridynamics models were added by Rezwanur Rahman and John
Foster (UTSA).

**Supporting info:**

* src/PERI: filenames -> commands
* `doc/PDF/PDLammps_overview.pdf <PDF/PDLammps_overview.pdf>`_
* `doc/PDF/PDLammps_EPS.pdf <PDF/PDLammps_EPS.pdf>`_
* `doc/PDF/PDLammps_VES.pdf <PDF/PDLammps_VES.pdf>`_
* :doc:`atom_style peri <atom_style>`
* :doc:`pair_style peri/\* <pair_peri>`
* :doc:`compute damage/atom <compute_damage_atom>`
* :doc:`compute plasticity/atom <compute_plasticity_atom>`
* examples/peri
* https://lammps.sandia.gov/movies.html#peri

----------

.. _PKG-PLUGIN:

PLUGIN package
--------------

**Contents:**

A :doc:`plugin <plugin>` command that can load and unload several
kind of styles in LAMMPS from shared object files at runtime without
having to recompile and relink LAMMPS.

**Authors:** Axel Kohlmeyer (Temple U)

**Supporting info:**

* src/PLUGIN: filenames -> commands
* :doc:`plugin command <plugin>`
* :doc:`Information on writing plugins <Developer_plugins>`
* examples/plugin

----------

.. _PKG-POEMS:

POEMS package
-------------

**Contents:**

A fix that wraps the Parallelizable Open source Efficient Multibody
Software (POEMS) library, which is able to simulate the dynamics of
articulated body systems.  These are systems with multiple rigid
bodies (collections of particles) whose motion is coupled by
connections at hinge points.

**Author:** Rudra Mukherjee (JPL) while at RPI.

**Install:**

This package has :ref:`specific installation instructions <poems>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/POEMS: filenames -> commands
* src/POEMS/README
* lib/poems/README
* :doc:`fix poems <fix_poems>`
* examples/rigid

----------

.. _PKG-PYTHON:

PYTHON package
--------------

**Contents:**

A :doc:`python <python>` command which allow you to execute Python code
from a LAMMPS input script.  The code can be in a separate file or
embedded in the input script itself.  See the :doc:`Python call <Python_call>` page for an overview of using Python from
LAMMPS in this manner and all the :doc:`Python <Python_head>` manual pages
for other ways to use LAMMPS and Python together.

.. note::

   Building with the PYTHON package assumes you have a Python
   shared library available on your system, which needs to be a Python 2
   version, 2.6 or later.  Python 3 is not yet supported.  See the
   lib/python/README for more details.

**Install:**

This package has :ref:`specific installation instructions <python>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/PYTHON: filenames -> commands
* :doc:`Python call <Python_head>`
* lib/python/README
* examples/python

----------

.. _PKG-QEQ:

QEQ package
-----------

**Contents:**

Several fixes for performing charge equilibration (QEq) via different
algorithms.  These can be used with pair styles that perform QEq as
part of their formulation.

**Supporting info:**

* src/QEQ: filenames -> commands
* :doc:`fix qeq/\* <fix_qeq>`
* examples/qeq
* examples/streitz

----------

.. _PKG-REPLICA:

REPLICA package
---------------

**Contents:**

A collection of multi-replica methods which can be used when running
multiple LAMMPS simulations (replicas).  See the :doc:`Howto replica <Howto_replica>` page for an overview of how to run
multi-replica simulations in LAMMPS.  Methods in the package include
nudged elastic band (NEB), parallel replica dynamics (PRD),
temperature accelerated dynamics (TAD), parallel tempering, and a
verlet/split algorithm for performing long-range Coulombics on one set
of processors, and the remainder of the force field calculation on
another set.

**Supporting info:**

* src/REPLICA: filenames -> commands
* :doc:`Howto replica <Howto_replica>`
* :doc:`neb <neb>`
* :doc:`prd <prd>`
* :doc:`tad <tad>`
* :doc:`temper <temper>`,
* :doc:`run_style verlet/split <run_style>`
* examples/neb
* examples/prd
* examples/tad

----------

.. _PKG-RIGID:

RIGID package
-------------

**Contents:**

Fixes which enforce rigid constraints on collections of atoms or
particles.  This includes SHAKE and RATTLE, as well as various
rigid-body integrators for a few large bodies or many small bodies.
Also several computes which calculate properties of rigid bodies.

**Supporting info:**

* src/RIGID: filenames -> commands
* :doc:`compute erotate/rigid <compute_erotate_rigid>`
* :doc:`fix shake <fix_shake>`
* :doc:`fix rattle <fix_shake>`
* :doc:`fix rigid/\* <fix_rigid>`
* examples/ASPHERE
* examples/rigid
* bench/in.rhodo
* https://lammps.sandia.gov/movies.html#box
* https://lammps.sandia.gov/movies.html#star

----------

.. _PKG-SHOCK:

SHOCK package
-------------

**Contents:**

Fixes for running impact simulations where a shock-wave passes through
a material.

**Supporting info:**

* src/SHOCK: filenames -> commands
* :doc:`fix append/atoms <fix_append_atoms>`
* :doc:`fix msst <fix_msst>`
* :doc:`fix nphug <fix_nphug>`
* :doc:`fix wall/piston <fix_wall_piston>`
* examples/hugoniostat
* examples/msst

----------

.. _PKG-SNAP:

SNAP package
------------

**Contents:**

A pair style for the spectral neighbor analysis potential (SNAP).
SNAP is methodology for deriving a highly accurate classical potential
fit to a large archive of quantum mechanical (DFT) data. Also several
computes which analyze attributes of the potential.

**Author:** Aidan Thompson (Sandia).

**Supporting info:**

* src/SNAP: filenames -> commands
* :doc:`pair_style snap <pair_snap>`
* :doc:`compute sna/atom <compute_sna_atom>`
* :doc:`compute snad/atom <compute_sna_atom>`
* :doc:`compute snav/atom <compute_sna_atom>`
* examples/snap

----------

.. _PKG-SPIN:

SPIN package
------------

**Contents:**

Model atomic magnetic spins classically, coupled to atoms moving in
the usual manner via MD.  Various pair, fix, and compute styles.

**Author:** Julien Tranchida (Sandia).

**Supporting info:**

* src/SPIN: filenames -> commands
* :doc:`Howto spins <Howto_spins>`
* :doc:`pair_style spin/dipole/cut <pair_spin_dipole>`
* :doc:`pair_style spin/dipole/long <pair_spin_dipole>`
* :doc:`pair_style spin/dmi <pair_spin_dmi>`
* :doc:`pair_style spin/exchange <pair_spin_exchange>`
* :doc:`pair_style spin/exchange/biquadratic <pair_spin_exchange>`
* :doc:`pair_style spin/magelec <pair_spin_magelec>`
* :doc:`pair_style spin/neel <pair_spin_neel>`
* :doc:`fix nve/spin <fix_nve_spin>`
* :doc:`fix langevin/spin <fix_langevin_spin>`
* :doc:`fix precession/spin <fix_precession_spin>`
* :doc:`compute spin <compute_spin>`
* :doc:`neb/spin <neb_spin>`
* examples/SPIN

----------

.. _PKG-SRD:

SRD package
-----------

**Contents:**

A pair of fixes which implement the Stochastic Rotation Dynamics (SRD)
method for coarse-graining of a solvent, typically around large
colloidal particles.

**Supporting info:**

* src/SRD: filenames -> commands
* :doc:`fix srd <fix_srd>`
* :doc:`fix wall/srd <fix_wall_srd>`
* examples/srd
* examples/ASPHERE
* https://lammps.sandia.gov/movies.html#tri
* https://lammps.sandia.gov/movies.html#line
* https://lammps.sandia.gov/movies.html#poly

----------

.. _PKG-VORONOI:

VORONOI package
---------------

**Contents:**

A compute command which calculates the Voronoi tesselation of a
collection of atoms by wrapping the `Voro++ library <voro-home_>`_.  This
can be used to calculate the local volume or each atoms or its near
neighbors.

.. _voro-home: http://math.lbl.gov/voro++

To use this package you must have the Voro++ library available on your
system.

**Author:** Daniel Schwen (INL) while at LANL.  The open-source Voro++
library was written by Chris Rycroft (Harvard U) while at UC Berkeley
and LBNL.

**Install:**

This package has :ref:`specific installation instructions <voronoi>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/VORONOI: filenames -> commands
* src/VORONOI/README
* lib/voronoi/README
* :doc:`compute voronoi/atom <compute_voronoi_atom>`
* examples/voronoi

----------

.. _PKG-USER-ADIOS:

USER-ADIOS package
------------------

**Contents:**

ADIOS is a high-performance I/O library. This package implements the
:doc:`dump atom/adios <dump_adios>`, :doc:`dump custom/adios <dump_adios>` and
:doc:`read_dump ... format adios <read_dump>`
commands to write and read data using the ADIOS library.

**Authors:** Norbert Podhorszki (ORNL) from the ADIOS developer team.

**Install:**

This package has :ref:`specific installation instructions <user-adios>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-ADIOS: filenames -> commands
* src/USER-ADIOS/README
* examples/USER/adios
* https://github.com/ornladios/ADIOS2
* :doc:`dump atom/adios <dump_adios>`
* :doc:`dump custom/adios <dump_adios>`
* :doc:`read_dump <read_dump>`

----------

.. _PKG-USER-ATC:

USER-ATC package
----------------

**Contents:**

ATC stands for atoms-to-continuum.  This package implements a :doc:`fix atc <fix_atc>` command to either couple molecular dynamics with
continuum finite element equations or perform on-the-fly conversion of
atomic information to continuum fields.

**Authors:** Reese Jones, Jeremy Templeton, Jon Zimmerman (Sandia).

**Install:**

This package has :ref:`specific installation instructions <user-atc>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-ATC: filenames -> commands
* src/USER-ATC/README
* :doc:`fix atc <fix_atc>`
* examples/USER/atc
* https://lammps.sandia.gov/pictures.html#atc

----------

.. _PKG-USER-AWPMD:

USER-AWPMD package
------------------

**Contents:**

AWPMD stands for Antisymmetrized Wave Packet Molecular Dynamics.  This
package implements an atom, pair, and fix style which allows electrons
to be treated as explicit particles in a classical molecular dynamics
model.

**Author:** Ilya Valuev (JIHT, Russia).

**Install:**

This package has :ref:`specific installation instructions <user-awpmd>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-AWPMD: filenames -> commands
* src/USER-AWPMD/README
* :doc:`pair_style awpmd/cut <pair_awpmd>`
* examples/USER/awpmd

----------

.. _PKG-USER-BOCS:

USER-BOCS package
-----------------

**Contents:**

This package provides :doc:`fix bocs <fix_bocs>`, a modified version
of :doc:`fix npt <fix_nh>` which includes the pressure correction to
the barostat as outlined in:

N. J. H. Dunn and W. G. Noid, "Bottom-up coarse-grained models that
accurately describe the structure, pressure, and compressibility of
molecular liquids," J. Chem. Phys. 143, 243148 (2015).

**Authors:** Nicholas J. H. Dunn and Michael R. DeLyser (The
Pennsylvania State University)

**Supporting info:**

The USER-BOCS user package for LAMMPS is part of the BOCS software package:
`https://github.com/noid-group/BOCS <https://github.com/noid-group/BOCS>`_

See the following reference for information about the entire package:

Dunn, NJH; Lebold, KM; DeLyser, MR; Rudzinski, JF; Noid, WG.
"BOCS: Bottom-Up Open-Source Coarse-Graining Software."
J. Phys. Chem. B. 122, 13, 3363-3377 (2018).

Example inputs are in the examples/USER/bocs folder.

----------

.. _PKG-USER-CGDNA:

USER-CGDNA package
------------------

**Contents:**

Several pair styles, bond styles, and integration fixes for coarse-grained
modelling of single- and double-stranded DNA and RNA based on the oxDNA and
oxRNA model of Doye, Louis and Ouldridge. The package includes Langevin-type
rigid-body integrators with improved stability.

**Author:** Oliver Henrich (University of Strathclyde, Glasgow).

**Supporting info:**

* src/USER-CGDNA: filenames -> commands
* /src/USER-CGDNA/README
* :doc:`pair_style oxdna/\* <pair_oxdna>`
* :doc:`pair_style oxdna2/\* <pair_oxdna2>`
* :doc:`pair_style oxrna2/\* <pair_oxrna2>`
* :doc:`bond_style oxdna/\* <bond_oxdna>`
* :doc:`bond_style oxdna2/\* <bond_oxdna>`
* :doc:`bond_style oxrna2/\* <bond_oxdna>`
* :doc:`fix nve/dotc/langevin <fix_nve_dotc_langevin>`

----------

.. _PKG-USER-CGSDK:

USER-CGSDK package
------------------

**Contents:**

Several pair styles and an angle style which implement the
coarse-grained SDK model of Shinoda, DeVane, and Klein which enables
simulation of ionic liquids, electrolytes, lipids and charged amino
acids.

**Author:** Axel Kohlmeyer (Temple U).

**Supporting info:**

* src/USER-CGSDK: filenames -> commands
* src/USER-CGSDK/README
* :doc:`pair_style lj/sdk/\* <pair_sdk>`
* :doc:`angle_style sdk <angle_sdk>`
* examples/USER/cgsdk
* https://lammps.sandia.gov/pictures.html#cg

----------

.. _PKG-USER-COLVARS:

USER-COLVARS package
--------------------

**Contents:**

COLVARS stands for collective variables, which can be used to
implement various enhanced sampling methods, including Adaptive
Biasing Force, Metadynamics, Steered MD, Umbrella Sampling and
Restraints.  A :doc:`fix colvars <fix_colvars>` command is implemented
which wraps a COLVARS library, which implements these methods.
simulations.

**Authors:** The COLVARS library is written and maintained by
Giacomo Fiorin (ICMS, Temple University, Philadelphia, PA, USA)
and Jerome Henin (LISM, CNRS, Marseille, France), originally for
the NAMD MD code, but with portability in mind.  Axel Kohlmeyer
(Temple U) provided the interface to LAMMPS.

**Install:**

This package has :ref:`specific installation instructions <user-colvars>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-COLVARS: filenames -> commands
* `doc/PDF/colvars-refman-lammps.pdf <PDF/colvars-refman-lammps.pdf>`_
* src/USER-COLVARS/README
* lib/colvars/README
* :doc:`fix colvars <fix_colvars>`
* examples/USER/colvars

----------

.. _PKG-USER-PACE:

USER-PACE package
-------------------

**Contents:**

A pair style for the Atomic Cluster Expansion potential (ACE).
ACE is a methodology for deriving a highly accurate classical potential
fit to a large archive of quantum mechanical (DFT) data. The USER-PACE
package provides an efficient implementation for running simulations
with ACE potentials.

**Authors:**

This package was written by Yury Lysogorskiy^1,
Cas van der Oord^2, Anton Bochkarev^1,
Sarath Menon^1, Matteo Rinaldi^1, Thomas Hammerschmidt^1, Matous Mrovec^1,
Aidan Thompson^3, Gabor Csanyi^2, Christoph Ortner^4, Ralf Drautz^1.

 ^1: Ruhr-University Bochum, Bochum, Germany

 ^2: University of Cambridge, Cambridge, United Kingdom

 ^3: Sandia National Laboratories, Albuquerque, New Mexico, USA

 ^4: University of British Columbia, Vancouver, BC, Canada

**Install:**

This package has :ref:`specific installation instructions <user-pace>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-PACE: filenames -> commands
* :doc:`pair_style pace <pair_pace>`
* examples/USER/pace

----------

.. _PKG-USER-PLUMED:

USER-PLUMED package
-------------------

**Contents:**

The fix plumed command allows you to use the PLUMED free energy plugin
for molecular dynamics to analyze and bias your LAMMPS trajectory on
the fly.  The PLUMED library is called from within the LAMMPS input
script by using the :doc:`fix plumed <fix_plumed>` command.

**Authors:** The :ref:`PLUMED library <PLUMED>` is written and maintained by
Massimilliano Bonomi, Giovanni Bussi, Carlo Camiloni and Gareth
Tribello.

.. _PLUMED: https://www.plumed.org

**Install:**

This package has :ref:`specific installation instructions <user-plumed>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-PLUMED/README
* lib/plumed/README
* :doc:`fix plumed <fix_plumed>`
* examples/USER/plumed

----------

.. _PKG-USER-DIFFRACTION:

USER-DIFFRACTION package
------------------------

**Contents:**

Two computes and a fix for calculating x-ray and electron diffraction
intensities based on kinematic diffraction theory.

**Author:** Shawn Coleman while at the U Arkansas.

**Supporting info:**

* src/USER-DIFFRACTION: filenames -> commands
* :doc:`compute saed <compute_saed>`
* :doc:`compute xrd <compute_xrd>`
* :doc:`fix saed/vtk <fix_saed_vtk>`
* examples/USER/diffraction

----------

.. _PKG-USER-DPD:

USER-DPD package
----------------

**Contents:**

DPD stands for dissipative particle dynamics.  This package implements
coarse-grained DPD-based models for energetic, reactive molecular
crystalline materials.  It includes many pair styles specific to these
systems, including for reactive DPD, where each particle has internal
state for multiple species and a coupled set of chemical reaction ODEs
are integrated each timestep.  Highly accurate time integrators for
isothermal, isoenergetic, isobaric and isenthalpic conditions are
included.  These enable long timesteps via the Shardlow splitting
algorithm.

**Authors:** Jim Larentzos (ARL), Tim Mattox (Engility Corp), and John
Brennan (ARL).

**Supporting info:**

* src/USER-DPD: filenames -> commands
* /src/USER-DPD/README
* :doc:`compute dpd <compute_dpd>`
* :doc:`compute dpd/atom <compute_dpd_atom>`
* :doc:`fix eos/cv <fix_eos_table>`
* :doc:`fix eos/table <fix_eos_table>`
* :doc:`fix eos/table/rx <fix_eos_table_rx>`
* :doc:`fix shardlow <fix_shardlow>`
* :doc:`fix rx <fix_rx>`
* :doc:`pair_style table/rx <pair_table_rx>`
* :doc:`pair_style dpd/fdt <pair_dpd_fdt>`
* :doc:`pair_style dpd/fdt/energy <pair_dpd_fdt>`
* :doc:`pair_style exp6/rx <pair_exp6_rx>`
* :doc:`pair_style multi/lucy <pair_multi_lucy>`
* :doc:`pair_style multi/lucy/rx <pair_multi_lucy_rx>`
* examples/USER/dpd

----------

.. _PKG-USER-DRUDE:

USER-DRUDE package
------------------

**Contents:**

Fixes, pair styles, and a compute to simulate thermalized Drude
oscillators as a model of polarization.  See the :doc:`Howto drude <Howto_drude>` and :doc:`Howto drude2 <Howto_drude2>` pages
for an overview of how to use the package.  There are auxiliary tools
for using this package in tools/drude.

**Authors:** Alain Dequidt (U Clermont Auvergne), Julien
Devemy (CNRS), and Agilio Padua (ENS de Lyon).

**Supporting info:**

* src/USER-DRUDE: filenames -> commands
* :doc:`Howto drude <Howto_drude>`
* :doc:`Howto drude2 <Howto_drude2>`
* :doc:`Howto polarizable <Howto_polarizable>`
* src/USER-DRUDE/README
* :doc:`fix drude <fix_drude>`
* :doc:`fix drude/transform/\* <fix_drude_transform>`
* :doc:`compute temp/drude <compute_temp_drude>`
* :doc:`pair_style thole <pair_thole>`
* :doc:`pair_style lj/cut/thole/long <pair_thole>`
* examples/USER/drude
* tools/drude

----------

.. _PKG-USER-EFF:

USER-EFF package
----------------

**Contents:**

EFF stands for electron force field which allows a classical MD code
to model electrons as particles of variable radius.  This package
contains atom, pair, fix and compute styles which implement the eFF as
described in A. Jaramillo-Botero, J. Su, Q. An, and W.A. Goddard III,
JCC, 2010.  The eFF potential was first introduced by Su and Goddard,
in 2007.  There are auxiliary tools for using this package in
tools/eff; see its README file.

**Author:** Andres Jaramillo-Botero (CalTech).

**Supporting info:**

* src/USER-EFF: filenames -> commands
* src/USER-EFF/README
* :doc:`atom_style electron <atom_style>`
* :doc:`fix nve/eff <fix_nve_eff>`
* :doc:`fix nvt/eff <fix_nh_eff>`
* :doc:`fix npt/eff <fix_nh_eff>`
* :doc:`fix langevin/eff <fix_langevin_eff>`
* :doc:`compute temp/eff <compute_temp_eff>`
* :doc:`pair_style eff/cut <pair_eff>`
* :doc:`pair_style eff/inline <pair_eff>`
* examples/USER/eff
* tools/eff/README
* tools/eff
* https://lammps.sandia.gov/movies.html#eff

----------

.. _PKG-USER-FEP:

USER-FEP package
----------------

**Contents:**

FEP stands for free energy perturbation.  This package provides
methods for performing FEP simulations by using a :doc:`fix adapt/fep <fix_adapt_fep>` command with soft-core pair potentials,
which have a "soft" in their style name.  There are auxiliary tools
for using this package in tools/fep; see its README file.

**Author:** Agilio Padua (ENS de Lyon)

**Supporting info:**

* src/USER-FEP: filenames -> commands
* src/USER-FEP/README
* :doc:`fix adapt/fep <fix_adapt_fep>`
* :doc:`compute fep <compute_fep>`
* :doc:`pair_style \*/soft <pair_fep_soft>`
* examples/USER/fep
* tools/fep/README
* tools/fep

----------

.. _PKG-USER-H5MD:

USER-H5MD package
-----------------

**Contents:**

H5MD stands for HDF5 for MD.  `HDF5 <HDF5_>`_ is a portable, binary,
self-describing file format, used by many scientific simulations.
H5MD is a format for molecular simulations, built on top of HDF5.
This package implements a :doc:`dump h5md <dump_h5md>` command to output
LAMMPS snapshots in this format.

.. _HDF5: http://www.hdfgroup.org/HDF5

To use this package you must have the HDF5 library available on your
system.

**Author:** Pierre de Buyl (KU Leuven) created both the package and the
H5MD format.

**Install:**

This package has :ref:`specific installation instructions <user-h5md>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-H5MD: filenames -> commands
* src/USER-H5MD/README
* lib/h5md/README
* :doc:`dump h5md <dump_h5md>`

----------

.. _PKG-USER-INTEL:

USER-INTEL package
------------------

**Contents:**

Dozens of pair, fix, bond, angle, dihedral, improper, and kspace
styles which are optimized for Intel CPUs and KNLs (Knights Landing).
All of them have an "intel" in their style name.  The
:doc:`USER-INTEL package <Speed_intel>` page gives details of what hardware and
compilers are required on your system, and how to build and use this
package.  Its styles can be invoked at run time via the "-sf intel" or
"-suffix intel" :doc:`command-line switches <Run_options>`.  Also see
the :ref:`KOKKOS <PKG-KOKKOS>`, :ref:`OPT <PKG-OPT>`, and :ref:`USER-OMP <PKG-USER-OMP>` packages,
which have styles optimized for CPUs and KNLs.

You need to have an Intel compiler, version 14 or higher to take full
advantage of this package. While compilation with GNU compilers is
supported, performance will be sub-optimal.

.. note::

   the USER-INTEL package contains styles that require using the
   -restrict flag, when compiling with Intel compilers.

**Author:** Mike Brown (Intel).

**Install:**

This package has :ref:`specific installation instructions <user-intel>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-INTEL: filenames -> commands
* src/USER-INTEL/README
* :doc:`Accelerator packages <Speed_packages>`
* :doc:`USER-INTEL package <Speed_intel>`
* :doc:`Section 2.6 -sf intel <Run_options>`
* :doc:`Section 2.6 -pk intel <Run_options>`
* :doc:`package intel <package>`
* Search the :doc:`commands <Commands_all>` pages (:doc:`fix <Commands_fix>`, :doc:`compute <Commands_compute>`,
  :doc:`pair <Commands_pair>`, :doc:`bond, angle, dihedral, improper <Commands_bond>`, :doc:`kspace <Commands_kspace>`) for styles followed by (i)
* src/USER-INTEL/TEST
* `Benchmarks page <https://lammps.sandia.gov/bench.html>`_ of web site

----------

.. _PKG-USER-LB:

USER-LB package
---------------

**Contents:**

Fixes which implement a background Lattice-Boltzmann (LB) fluid, which
can be used to model MD particles influenced by hydrodynamic forces.

**Authors:** Frances Mackay and Colin Denniston (University of Western
Ontario).

**Supporting info:**

* src/USER-LB: filenames -> commands
* src/USER-LB/README
* :doc:`fix lb/fluid <fix_lb_fluid>`
* :doc:`fix lb/momentum <fix_lb_momentum>`
* :doc:`fix lb/viscous <fix_lb_viscous>`
* examples/USER/lb

----------

.. _PKG-USER-MGPT:

USER-MGPT package
-----------------

**Contents:**

A pair style which provides a fast implementation of the quantum-based
MGPT multi-ion potentials.  The MGPT or model GPT method derives from
first-principles DFT-based generalized pseudopotential theory (GPT)
through a series of systematic approximations valid for mid-period
transition metals with nearly half-filled d bands.  The MGPT method
was originally developed by John Moriarty at LLNL.  The pair style in
this package calculates forces and energies using an optimized
matrix-MGPT algorithm due to Tomas Oppelstrup at LLNL.

**Authors:** Tomas Oppelstrup and John Moriarty (LLNL).

**Supporting info:**

* src/USER-MGPT: filenames -> commands
* src/USER-MGPT/README
* :doc:`pair_style mgpt <pair_mgpt>`
* examples/USER/mgpt

----------

.. _PKG-USER-MISC:

USER-MISC package
-----------------

**Contents:**

A potpourri of (mostly) unrelated features contributed to LAMMPS by
users.  Each feature is a single fix, compute, pair, bond, angle,
dihedral, improper, or command style.

**Authors:** The author for each style in the package is listed in the
src/USER-MISC/README file.

**Supporting info:**

* src/USER-MISC: filenames -> commands
* src/USER-MISC/README
* one page per individual command listed in src/USER-MISC/README
* examples/USER/misc

----------

.. _PKG-USER-MANIFOLD:

USER-MANIFOLD package
---------------------

**Contents:**

Several fixes and a "manifold" class which enable simulations of
particles constrained to a manifold (a 2D surface within the 3D
simulation box).  This is done by applying the RATTLE constraint
algorithm to formulate single-particle constraint functions
g(xi,yi,zi) = 0 and their derivative (i.e. the normal of the manifold)
n = grad(g).

**Author:** Stefan Paquay (until 2017: Eindhoven University of
Technology (TU/e), The Netherlands; since 2017: Brandeis University,
Waltham, MA, USA)

**Supporting info:**

* src/USER-MANIFOLD: filenames -> commands
* src/USER-MANIFOLD/README
* :doc:`Howto manifold <Howto_manifold>`
* :doc:`fix manifoldforce <fix_manifoldforce>`
* :doc:`fix nve/manifold/rattle <fix_nve_manifold_rattle>`
* :doc:`fix nvt/manifold/rattle <fix_nvt_manifold_rattle>`
* examples/USER/manifold
* https://lammps.sandia.gov/movies.html#manifold

----------

.. _PKG-USER-MEAMC:

USER-MEAMC package
------------------

**Contents:**

A pair style for the modified embedded atom (MEAM) potential
translated from the Fortran version in the (obsolete) MEAM package
to plain C++. The USER-MEAMC fully replaces the MEAM package, which
has been removed from LAMMPS after the 12 December 2018 version.

**Author:** Sebastian Huetter, (Otto-von-Guericke University Magdeburg)
based on the Fortran version of Greg Wagner (Northwestern U) while at
Sandia.

**Supporting info:**

* src/USER-MEAMC: filenames -> commands
* src/USER-MEAMC/README
* :doc:`pair_style meam/c <pair_meamc>`
* examples/meamc

----------

.. _PKG-USER-MESODPD:

USER-MESODPD package
--------------------

**Contents:**

Several extensions of the dissipative particle dynamics (DPD)
method.  Specifically, energy-conserving DPD (eDPD) that can model
non-isothermal processes, many-body DPD (mDPD) for simulating
vapor-liquid coexistence, and transport DPD (tDPD) for modeling
advection-diffusion-reaction systems. The equations of motion of these
DPD extensions are integrated through a modified velocity-Verlet (MVV)
algorithm.

**Author:** Zhen Li (Division of Applied Mathematics, Brown University)

**Supporting info:**

* src/USER-MESODPD: filenames -> commands
* src/USER-MESODPD/README
* :doc:`atom_style edpd <atom_style>`
* :doc:`pair_style edpd <pair_mesodpd>`
* :doc:`pair_style mdpd <pair_mesodpd>`
* :doc:`pair_style tdpd <pair_mesodpd>`
* :doc:`fix mvv/dpd <fix_mvv_dpd>`
* examples/USER/mesodpd
* https://lammps.sandia.gov/movies.html#mesodpd

* examples/USER/meso
* http://lammps.sandia.gov/movies.html#mesodpd

----------

.. _PKG-USER-MESONT:

USER-MESONT package
-------------------

**Contents:**

USER-MESONT is a LAMMPS package for simulation of nanomechanics of
nanotubes (NTs). The model is based on a coarse-grained representation
of NTs as "flexible cylinders" consisting of a variable number of
segments. Internal interactions within a NT and the van der Waals
interaction between the tubes are described by a mesoscopic force field
designed and parameterized based on the results of atomic-level
molecular dynamics simulations. The description of the force field is
provided in the papers listed below. This package contains two
independent implementations of this model: :doc:`pair_style mesocnt
<pair_mesocnt>` is a (minimal) C++ implementation, and :doc:`pair_style
mesont/tpm <pair_mesont_tpm>` is a more general and feature rich
implementation based on a Fortran library in the ``lib/mesont`` folder.

**Download of potential files:**

The potential files for these pair styles are *very* large and thus
are not included in the regular downloaded packages of LAMMPS or the
git repositories.  Instead, they will be automatically downloaded
from a web server when the package is installed for the first time.

**Authors of the *mesont* styles:**

Maxim V. Shugaev (University of Virginia), Alexey N. Volkov (University of Alabama), Leonid V. Zhigilei (University of Virginia)

**Author of the *mesocnt* pair style:**
Philipp Kloza (U Cambridge)

**Supporting info:**

* src/USER-MESONT: filenames -> commands
* src/USER-MESONT/README
* :doc:`atom_style mesont <atom_style>`
* :doc:`pair_style mesont/tpm <pair_mesont_tpm>`
* :doc:`compute mesont <compute_mesont>`
* :doc:`pair_style mesocnt <pair_mesocnt>`
* examples/USER/mesont
* tools/mesont

----------

.. _PKG-USER-MOFFF:

USER-MOFFF package
------------------

**Contents:**

Pair, angle and improper styles needed to employ the MOF-FF
force field by Schmid and coworkers with LAMMPS.
MOF-FF is a first principles derived force field with the primary aim
to simulate MOFs and related porous framework materials, using spherical
Gaussian charges. It is described in S. Bureekaew et al., Phys. Stat. Sol. B
2013, 250, 1128-1141.
For the usage of MOF-FF see the example in the example directory as
well as the `MOF+ <MOFplus_>`_ website.

.. _MOFplus: https://www.mofplus.org/content/show/MOF-FF

**Author:** Hendrik Heenen (Technical U of Munich),
Rochus Schmid (Ruhr-University Bochum).

**Supporting info:**

* src/USER-MOFFF: filenames -> commands
* src/USER-MOFFF/README
* :doc:`pair_style buck6d/coul/gauss <pair_buck6d_coul_gauss>`
* :doc:`angle_style class2 <angle_class2>`
* :doc:`angle_style cosine/buck6d <angle_cosine_buck6d>`
* :doc:`improper_style inversion/harmonic <improper_inversion_harmonic>`
* examples/USER/mofff

----------

.. _PKG-USER-MOLFILE:

USER-MOLFILE package
--------------------

**Contents:**

A :doc:`dump molfile <dump_molfile>` command which uses molfile plugins
that are bundled with the `VMD <vmd-home_>`_
molecular visualization and analysis program, to enable LAMMPS to dump
snapshots in formats compatible with various molecular simulation
tools.

To use this package you must have the desired VMD plugins available on
your system.

Note that this package only provides the interface code, not the
plugins themselves, which will be accessed when requesting a specific
plugin via the :doc:`dump molfile <dump_molfile>` command.  Plugins can
be obtained from a VMD installation which has to match the platform
that you are using to compile LAMMPS for. By adding plugins to VMD,
support for new file formats can be added to LAMMPS (or VMD or other
programs that use them) without having to re-compile the application
itself.  More information about the VMD molfile plugins can be found
at
`http://www.ks.uiuc.edu/Research/vmd/plugins/molfile <http://www.ks.uiuc.edu/Research/vmd/plugins/molfile>`_.

**Author:** Axel Kohlmeyer (Temple U).

**Install:**

This package has :ref:`specific installation instructions <user-molfile>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-MOLFILE: filenames -> commands
* src/USER-MOLFILE/README
* lib/molfile/README
* :doc:`dump molfile <dump_molfile>`

----------

.. _PKG-USER-NETCDF:

USER-NETCDF package
-------------------

**Contents:**

Dump styles for writing NetCDF formatted dump files.  NetCDF is a
portable, binary, self-describing file format developed on top of
HDF5. The file contents follow the AMBER NetCDF trajectory conventions
(http://ambermd.org/netcdf/nctraj.xhtml), but include extensions.

To use this package you must have the NetCDF library available on your
system.

Note that NetCDF files can be directly visualized with the following
tools:

* `Ovito <ovito_>`_ (Ovito supports the AMBER convention and the extensions mentioned above)
* `VMD <vmd-home_>`_

.. _ovito: http://www.ovito.org

.. _vmd-home: https://www.ks.uiuc.edu/Research/vmd/

**Author:** Lars Pastewka (Karlsruhe Institute of Technology).

**Install:**

This package has :ref:`specific installation instructions <user-netcdf>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-NETCDF: filenames -> commands
* src/USER-NETCDF/README
* lib/netcdf/README
* :doc:`dump netcdf <dump_netcdf>`

----------

.. _PKG-USER-OMP:

USER-OMP package
----------------

**Contents:**

Hundreds of pair, fix, compute, bond, angle, dihedral, improper, and
kspace styles which are altered to enable threading on many-core CPUs
via OpenMP directives.  All of them have an "omp" in their style name.
The :doc:`USER-OMP package <Speed_omp>` page gives details of what hardware
and compilers are required on your system, and how to build and use
this package.  Its styles can be invoked at run time via the "-sf omp"
or "-suffix omp" :doc:`command-line switches <Run_options>`.  Also see
the :ref:`KOKKOS <PKG-KOKKOS>`, :ref:`OPT <PKG-OPT>`, and :ref:`USER-INTEL <PKG-USER-INTEL>`
packages, which have styles optimized for CPUs.

**Author:** Axel Kohlmeyer (Temple U).

.. note::

   To enable multi-threading support the compile flag "-fopenmp"
   and the link flag "-fopenmp" (for GNU compilers, you have to look up
   the equivalent flags for other compilers) must be used to build LAMMPS.
   When using Intel compilers, also the "-restrict" flag is required.
   The USER-OMP package can be compiled without enabling OpenMP; then
   all code will be compiled as serial and the only improvement over the
   regular styles are some data access optimization. These flags should
   be added to the CCFLAGS and LINKFLAGS lines of your Makefile.machine.
   See src/MAKE/OPTIONS/Makefile.omp for an example.

Once you have an appropriate Makefile.machine, you can
install/un-install the package and build LAMMPS in the usual manner:

**Install:**

This package has :ref:`specific installation instructions <user-omp>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-OMP: filenames -> commands
* src/USER-OMP/README
* :doc:`Accelerator packages <Speed_packages>`
* :doc:`USER-OMP package <Speed_omp>`
* :doc:`Section 2.6 -sf omp <Run_options>`
* :doc:`Section 2.6 -pk omp <Run_options>`
* :doc:`package omp <package>`
* Search the :doc:`commands <Commands_all>` pages (:doc:`fix <Commands_fix>`, :doc:`compute <Commands_compute>`,
  :doc:`pair <Commands_pair>`, :doc:`bond, angle, dihedral, improper <Commands_bond>`,
  :doc:`kspace <Commands_kspace>`) for styles followed by (o)
* `Benchmarks page <https://lammps.sandia.gov/bench.html>`_ of web site

----------

.. _PKG-USER-PHONON:

USER-PHONON package
-------------------

**Contents:**

A :doc:`fix phonon <fix_phonon>` command that calculates dynamical
matrices, which can then be used to compute phonon dispersion
relations, directly from molecular dynamics simulations.
And a :doc:`dynamical_matrix <dynamical_matrix>` as well as a
:doc:`third_order <third_order>` command to compute the dynamical matrix
and third order tensor from finite differences.

**Authors:** Ling-Ti Kong (Shanghai Jiao Tong University) for "fix phonon"
and Charlie Sievers (UC Davis) for "dynamical_matrix" and "third_order"

**Supporting info:**

* src/USER-PHONON: filenames -> commands
* src/USER-PHONON/README
* :doc:`fix phonon <fix_phonon>`
* :doc:`dynamical_matrix <dynamical_matrix>`
* :doc:`third_order <third_order>`
* examples/USER/phonon

----------

.. _PKG-USER-PTM:

USER-PTM package
----------------

**Contents:**

A :doc:`compute ptm/atom <compute_ptm_atom>` command that calculates
local structure characterization using the Polyhedral Template
Matching methodology.

**Author:** Peter Mahler Larsen (MIT).

**Supporting info:**

* src/USER-PTM: filenames not starting with ptm\_ -> commands
* src/USER-PTM: filenames starting with ptm\_ -> supporting code
* src/USER-PTM/LICENSE
* :doc:`compute ptm/atom <compute_ptm_atom>`

----------

.. _PKG-USER-QMMM:

USER-QMMM package
-----------------

**Contents:**

A :doc:`fix qmmm <fix_qmmm>` command which allows LAMMPS to be used as
the MM code in a QM/MM simulation.  This is currently only available
in combination with the `Quantum ESPRESSO <espresso_>`_ package.

.. _espresso: http://www.quantum-espresso.org

To use this package you must have Quantum ESPRESSO (QE) available on
your system and include its coupling library in the compilation and
then compile LAMMPS as a library.  For QM/MM calculations you then
build a custom binary with MPI support, that sets up 3 partitions with
MPI sub-communicators (for inter- and intra-partition communication)
and then calls the corresponding library interfaces on each partition
(2x LAMMPS and 1x QE).

The current implementation supports an ONIOM style mechanical coupling
and a multi-pole based electrostatic coupling to the Quantum ESPRESSO
plane wave DFT package.  The QM/MM interface has been written in a
manner that coupling to other QM codes should be possible without
changes to LAMMPS itself.

**Authors:** Axel Kohlmeyer (Temple U). Mariella Ippolito and Carlo Cavazzoni (CINECA, Italy)

**Install:**

This package has :ref:`specific installation instructions <user-qmmm>`
on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-QMMM: filenames -> commands
* src/USER-QMMM/README
* lib/qmmm/README
* :doc:`fix phonon <fix_phonon>`
* lib/qmmm/example-ec/README
* lib/qmmm/example-mc/README

----------

.. _PKG-USER-QTB:

USER-QTB package
----------------

**Contents:**

Two fixes which provide a self-consistent quantum treatment of
vibrational modes in a classical molecular dynamics simulation.  By
coupling the MD simulation to a colored thermostat, it introduces zero
point energy into the system, altering the energy power spectrum and
the heat capacity to account for their quantum nature. This is useful
when modeling systems at temperatures lower than their classical
limits or when temperatures ramp across the classical limits in a
simulation.

**Author:** Yuan Shen (Stanford U).

**Supporting info:**

* src/USER-QTB: filenames -> commands
* src/USER-QTB/README
* :doc:`fix qtb <fix_qtb>`
* :doc:`fix qbmsst <fix_qbmsst>`
* examples/USER/qtb

----------

.. _PKG-USER-QUIP:

USER-QUIP package
-----------------

**Contents:**

A :doc:`pair_style quip <pair_quip>` command which wraps the `QUIP libAtoms library <quip_>`_, which includes a variety of interatomic
potentials, including Gaussian Approximation Potential (GAP) models
developed by the Cambridge University group.

.. _quip: https://github.com/libAtoms/QUIP

To use this package you must have the QUIP libAtoms library available
on your system.

**Author:** Albert Bartok (Cambridge University)

**Install:**

This package has :ref:`specific installation instructions <user-quip>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-QUIP: filenames -> commands
* src/USER-QUIP/README
* :doc:`pair_style quip <pair_quip>`
* examples/USER/quip

----------

.. _PKG-USER-REACTION:

USER-REACTION package
---------------------

**Contents:**

This package allows for complex bond topology changes (reactions)
during a running MD simulation, when using classical force fields.
Topology changes are defined in pre- and post-reaction molecule
templates and can include creation and deletion of bonds, angles,
dihedrals, impropers, atom types, bond types, angle types, dihedral
types, improper types, and/or atomic charges. Other options currently
available include reaction constraints (e.g. angle and Arrhenius
constraints), deletion of reaction byproducts or other small
molecules, and chiral-sensitive reactions.

**Author:** Jacob R. Gissinger (CU Boulder) while at NASA Langley Research Center.

**Supporting info:**

* src/USER-REACTION: filenames -> commands
* src/USER-REACTION/README
* :doc:`fix bond/react <fix_bond_react>`
* examples/USER/reaction
* `2017 LAMMPS Workshop <https://lammps.sandia.gov/workshops/Aug17/pdf/gissinger.pdf>`_
* `2019 LAMMPS Workshop <https://lammps.sandia.gov/workshops/Aug19/talk_gissinger.pdf>`_
* reacter.org

----------

.. _PKG-USER-REAXC:

USER-REAXC package
------------------

**Contents:**

A pair style which implements the ReaxFF potential in C/C++.  ReaxFF
is a universal reactive force field.  See the src/USER-REAXC/README file
for more info on differences between the two packages.  Also two fixes
for monitoring molecules as bonds are created and destroyed.

**Author:** Hasan Metin Aktulga (MSU) while at Purdue University.

**Supporting info:**

* src/USER-REAXC: filenames -> commands
* src/USER-REAXC/README
* :doc:`pair_style reax/c <pair_reaxc>`
* :doc:`fix reax/c/bonds <fix_reaxc_bonds>`
* :doc:`fix reax/c/species <fix_reaxc_species>`
* examples/reax

----------

.. _PKG-USER-SCAFACOS:

USER-SCAFACOS package
---------------------

**Contents:**

A KSpace style which wraps the `ScaFaCoS Coulomb solver library <http://www.scafacos.de>`_ to compute long-range Coulombic
interactions.

To use this package you must have the ScaFaCoS library available on
your system.

**Author:** Rene Halver (JSC) wrote the scafacos LAMMPS command.

ScaFaCoS itself was developed by a consortium of German research
facilities with a BMBF (German Ministry of Science and Education)
funded project in 2009-2012. Participants of the consortium were the
Universities of Bonn, Chemnitz, Stuttgart, and Wuppertal as well as
the Forschungszentrum Juelich.

**Install:**

This package has :ref:`specific installation instructions <user-scafacos>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-SCAFACOS: filenames -> commands
* src/USER-SCAFACOS/README
* :doc:`kspace_style scafacos <kspace_style>`
* :doc:`kspace_modify <kspace_modify>`
* examples/USER/scafacos

----------

.. _PKG-USER-SDPD:

USER-SDPD package
-----------------

**Contents:**

A pair style for smoothed dissipative particle dynamics (SDPD), which
is an extension of smoothed particle hydrodynamics (SPH) to mesoscale
where thermal fluctuations are important (see the
:ref:`USER-SPH package <PKG-USER-SPH>`).
Also two fixes for moving and rigid body integration of SPH/SDPD particles
(particles of atom_style meso).

**Author:** Morteza Jalalvand (Institute for Advanced Studies in Basic
Sciences, Iran).

**Supporting info:**

* src/USER-SDPD: filenames -> commands
* src/USER-SDPD/README
* :doc:`pair_style sdpd/taitwater/isothermal <pair_sdpd_taitwater_isothermal>`
* :doc:`fix meso/move <fix_meso_move>`
* :doc:`fix rigid/meso <fix_rigid_meso>`
* examples/USER/sdpd

----------

.. _PKG-USER-SMD:

USER-SMD package
----------------

**Contents:**

An atom style, fixes, computes, and several pair styles which
implements smoothed Mach dynamics (SMD) for solids, which is a model
related to smoothed particle hydrodynamics (SPH) for liquids (see the
:ref:`USER-SPH package <PKG-USER-SPH>`).

This package solves solids mechanics problems via a state of the art
stabilized meshless method with hourglass control.  It can specify
hydrostatic interactions independently from material strength models,
i.e. pressure and deviatoric stresses are separated.  It provides many
material models (Johnson-Cook, plasticity with hardening,
Mie-Grueneisen, Polynomial EOS) and allows new material models to be
added.  It implements rigid boundary conditions (walls) which can be
specified as surface geometries from \*.STL files.

**Author:** Georg Ganzenmuller (Fraunhofer-Institute for High-Speed
Dynamics, Ernst Mach Institute, Germany).

**Install:**

This package has :ref:`specific installation instructions <user-smd>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-SMD: filenames -> commands
* src/USER-SMD/README
* doc/PDF/SMD_LAMMPS_userguide.pdf
* examples/USER/smd
* https://lammps.sandia.gov/movies.html#smd

----------

.. _PKG-USER-SMTBQ:

USER-SMTBQ package
------------------

**Contents:**

A pair style which implements a Second Moment Tight Binding model with
QEq charge equilibration (SMTBQ) potential for the description of
ionocovalent bonds in oxides.

**Authors:** Nicolas Salles, Emile Maras, Olivier Politano, and Robert
Tetot (LAAS-CNRS, France).

**Supporting info:**

* src/USER-SMTBQ: filenames -> commands
* src/USER-SMTBQ/README
* :doc:`pair_style smtbq <pair_smtbq>`
* examples/USER/smtbq

----------

.. _PKG-USER-SPH:

USER-SPH package
----------------

**Contents:**

An atom style, fixes, computes, and several pair styles which
implements smoothed particle hydrodynamics (SPH) for liquids.  See the
related :ref:`USER-SMD package <PKG-USER-SMD>` package for smooth Mach dynamics
(SMD) for solids.

This package contains ideal gas, Lennard-Jones equation of states,
Tait, and full support for complete (i.e. internal-energy dependent)
equations of state.  It allows for plain or Monaghans XSPH integration
of the equations of motion.  It has options for density continuity or
density summation to propagate the density field.  It has
:doc:`set <set>` command options to set the internal energy and density
of particles from the input script and allows the same quantities to
be output with thermodynamic output or to dump files via the :doc:`compute property/atom <compute_property_atom>` command.

**Author:** Georg Ganzenmuller (Fraunhofer-Institute for High-Speed
Dynamics, Ernst Mach Institute, Germany).

**Supporting info:**

* src/USER-SPH: filenames -> commands
* src/USER-SPH/README
* doc/PDF/SPH_LAMMPS_userguide.pdf
* examples/USER/sph
* https://lammps.sandia.gov/movies.html#sph

----------

.. _PKG-USER-TALLY:

USER-TALLY package
------------------

**Contents:**

Several compute styles that can be called when pairwise interactions
are calculated to tally information (forces, heat flux, energy,
stress, etc) about individual interactions.

**Author:** Axel Kohlmeyer (Temple U).

**Supporting info:**

* src/USER-TALLY: filenames -> commands
* src/USER-TALLY/README
* :doc:`compute \*/tally <compute_tally>`
* examples/USER/tally

----------

.. _PKG-USER-UEF:

USER-UEF package
----------------

**Contents:**

A fix style for the integration of the equations of motion under
extensional flow with proper boundary conditions, as well as several
supporting compute styles and an output option.

**Author:** David Nicholson (MIT).

**Supporting info:**

* src/USER-UEF: filenames -> commands
* src/USER-UEF/README
* :doc:`fix nvt/uef <fix_nh_uef>`
* :doc:`fix npt/uef <fix_nh_uef>`
* :doc:`compute pressure/uef <compute_pressure_uef>`
* :doc:`compute temp/uef <compute_temp_uef>`
* :doc:`dump cfg/uef <dump_cfg_uef>`
* examples/uef

----------

.. _PKG-USER-VTK:

USER-VTK package
----------------

**Contents:**

A :doc:`dump vtk <dump_vtk>` command which outputs snapshot info in the
`VTK format <vtk_>`_, enabling visualization by `Paraview <paraview_>`_ or
other visualization packages.

.. _vtk: http://www.vtk.org

.. _paraview: http://www.paraview.org

To use this package you must have VTK library available on your
system.

**Authors:** Richard Berger (JKU) and Daniel Queteschiner (DCS Computing).

**Install:**

This package has :ref:`specific installation instructions <user-vtk>` on the :doc:`Build extras <Build_extras>` page.

**Supporting info:**

* src/USER-VTK: filenames -> commands
* src/USER-VTK/README
* lib/vtk/README
* :doc:`dump vtk <dump_vtk>`

----------

.. _PKG-USER-YAFF:

USER-YAFF package
-----------------

**Contents:**

Some potentials that are also implemented in the Yet Another Force Field (`YAFF <yaff_>`_) code.
The expressions and their use are discussed in the following papers

* Vanduyfhuys et al., J. Comput. Chem., 36 (13), 1015-1027 (2015) `link <vanduyfhuys2015_>`_
* Vanduyfhuys et al., J. Comput. Chem., 39 (16), 999-1011 (2018) `link <vanduyfhuys2018_>`_

which discuss the `QuickFF <quickff_>`_ methodology.

.. _vanduyfhuys2015: https://doi.org/10.1002/jcc.23877
.. _vanduyfhuys2018: https://doi.org/10.1002/jcc.25173
.. _quickff: http://molmod.github.io/QuickFF
.. _yaff: https://github.com/molmod/yaff

**Author:** Steven Vandenbrande.

**Supporting info:**

* src/USER-YAFF/README
* :doc:`angle_style cross <angle_cross>`
* :doc:`angle_style mm3 <angle_mm3>`
* :doc:`bond_style mm3 <bond_mm3>`
* :doc:`improper_style distharm <improper_distharm>`
* :doc:`improper_style sqdistharm <improper_sqdistharm>`
* :doc:`pair_style mm3/switch3/coulgauss/long <pair_lj_switch3_coulgauss_long>`
* :doc:`pair_style lj/switch3/coulgauss/long <pair_lj_switch3_coulgauss_long>`
* examples/USER/yaff
