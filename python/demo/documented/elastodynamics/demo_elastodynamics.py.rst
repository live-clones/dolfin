
.. _demo_elastodynamics:

============================================
Time-integration of elastodynamics equation
============================================

This demo is implemented in a single Python file,
:download:`demo_elastodynamics.py`, which contains both the
variational forms and the solver.

This demo shows how to perform time integration of transient elastodynamics using the Newmark scheme.
In particular it demonstrates how to:

* formulate mass and damping forms of the elastodynamics equation
* obtain a mass matrix lumped at nodes
* implement the Newmark scheme and its influence on the solution
* perform efficient computation of stresses using ``LocalSolver``


----------------------------------------- 
Introduction and elastodynamics equation
-----------------------------------------

The elastodynamics equation combine the balance of linear momentum:

.. math::
   \nabla \cdot \sigma + \rho b = \rho \ddot{u}

where :math:`u` is the displacement vector field, :math:`\rho` the material density, 
:math:`b` a given body force and :math:`\sigma` the stress tensor which is related 
to the displacement through a constitutive equation. In the case of isotropic linearized elasticity, one has:

.. math::
   \sigma =\lambda \text{tr}(\varepsilon)\mathbb{1} + 2\mu\varepsilon

where :math:`\varepsilon = (\nabla u + (\nabla u)^T)/2` is the linearized strain tensor, :math:`\mathbb{1}` is the 
identity of second-rank tensors and :math:`\lambda=\dfrac{E\nu}{(1+\nu)(1-2\nu)},\mu=\dfrac{E}{2(1+\nu)}` are the 
Lamé coefficients given as functions of the Young modulus :math:`E` and the Poisson ratio :math:`\nu`.

The weak form is readily obtained by integrating by part the balance equation using a test function :math:`v\in V`
with :math:`V` being a suitable function space that satisfies the displacement boundary conditions:

.. math::
   \int_{\Omega} \rho \ddot{u}\cdot v \, {\rm d} x + \int_{\Omega} \sigma(u):\varepsilon(v) \, {\rm d} x = 
   \int_{\Omega} \rho b \cdot v  \, {\rm d} x + \int_{\partial\Omega} (\sigma\cdot n) \cdot v \, {\rm d} s \quad \text{for all } v\in V

The previous equation can be written as follows:

.. math::
   \text{Find }u\in V\text{ such that } m(\ddot{u},v) + k(u,v) = L(v) \quad \text{for all } v\in V

where :math:`m` is the symmetric bilinear form associated with the mass matrix and :math:`k` the one associated with the stiffness matrix.

After introducing the finite element space interpolation, one obtains the corresponding discretized evolution equation:

.. math::
   \text{Find }\{u\}\in\mathbb{R}^n\text{ such that } \{v\}^T[M]\{\ddot{u}\} + \{v\}^T[K]\{u\} = \{v\}^T\{F\} \quad \text{for all } \{v\}\in\mathbb{R}^n

which is a generalized :math:`n`-dof harmonic oscillator equation.

Quite often in structural dynamics, structures do not oscillate perfectly but lose energy through various dissipative mechanisms (friction with air or supports,
internal dissipation through plasticity, damage, etc.). Dissipative terms can be introduced at the level of the constitutive equation if these mechanisms are well
known but quite often it is not the case. Dissipation can then be modeled by adding an ad hoc damping term depending on the structure velocity :math:`\dot{u}` 
to the previous evolution equation:

.. math::
   \text{Find }u\in V\text{ such that } m(\ddot{u},v) + c(\dot{u},v) + k(u,v) = L(v) \quad \text{for all } v\in V

The damping form will be considered here as bilinear and symmetric, being therefore associated with a damping matrix :math:`[C]`.

-----------------------------------------
Time discretization using Newmark scheme
-----------------------------------------


---------------
Implementation
---------------

When little is known about the origin of damping in the structure, a popular choice for the damping matrix, known as *Rayleigh damping*, consists in using
a linear combination of the mass and stiffness matrix :math:`[C] = \alpha_M[M]+\alpha_K[K]` with two positive parameters :math:`\alpha_M,\alpha_K` which 
can be fitted against experimental measures for instance (usually by measuring the damping ratio of two natural modes of vibration).


After importing the relevant modules, the geometry of a beam of length :math:`L=20`
and rectangular section of size :math:`B\times H` with :math:`B=0.5, H=1` is first defined::

 from fenics import *
 import numpy as np

 L, B, H = 20., 0.5, 1.

 Nx = 200
 Ny = int(B/L*Nx)+1
 Nz = int(H/L*Nx)+1

 mesh = BoxMesh(Point(0.,0.,0.),Point(L,B,H), Nx, Ny, Nz)


Material parameters and elastic constitutive relations are classical (here we
take :math:`\nu=0`) and we also introduce the material density :math:`\rho` for
later definition of the mass matrix::

 E, nu = 1e5, 0.
 rho = 1e-3

 # Lame coefficient for constitutive relation
 mu = E/2./(1+nu)
 lmbda = E*nu/(1+nu)/(1-2*nu)

 def eps(v):
     return sym(grad(v))
 def sigma(v):
     dim = v.geometric_dimension()
     return 2.0*mu*eps(v) + lmbda*tr(eps(v))*Identity(dim)

Standard FunctionSpace is defined and boundary conditions correspond to a
fully clamped support at :math:`x=0`::

 V = VectorFunctionSpace(mesh, 'Lagrange', degree=1)
 u_ = TrialFunction(V)
 du = TestFunction(V)


 def left(x, on_boundary):
     return near(x[0],0.)

 bc = DirichletBC(V, Constant((0.,0.,0.)), left)


The system stiffness matrix :math:`[K]` and mass matrix :math:`[M]` are
respectively obtained from assembling the corresponding variational forms::

 k_form = inner(sigma(du),eps(u_))*dx
 l_form = Constant(1.)*u_[0]*dx
 K = PETScMatrix()
 b = PETScVector()
 assemble_system(k_form, l_form, bc, A_tensor=K, b_tensor=b)

 m_form = rho*dot(du,u_)*dx
 M = PETScMatrix()
 assemble(m_form, tensor=M)

Matrices :math:`[K]` and :math:`[M]` are first defined as PETSc Matrix and
forms are assembled into it to ensure that they have the right type.
Note that boundary conditions have been applied to the stiffness matrix using
``assemble_system`` so as to preserve symmetry (a dummy ``l_form`` and right-hand side
vector have been introduced to call this function).


Modal dynamic analysis consists in solving the following generalized
eigenvalue problem :math:`[K]\{U\}=\lambda[M]\{U\}` where the eigenvalue
is related to the eigenfrequency :math:`\lambda=\omega^2`. This problem
can be solved using the ``SLEPcEigenSolver``. ::

 eigensolver = SLEPcEigenSolver(K, M)
 eigensolver.parameters['problem_type'] = 'gen_hermitian'
 eigensolver.parameters["spectrum"] = "smallest real"
 eigensolver.parameters['spectral_transform'] = 'shift-and-invert'
 eigensolver.parameters['spectral_shift'] = 0.

The problem type is specified to be a generalized eigenvalue problem with
Hermitian matrices. By default, SLEPc computes the largest eigenvalues, here
we instead look for the smallest eigenvalues (they should all be real). To
improve convergence of the eigensolver for finding the smallest eigenvalues
(by default it computes the largest ones), a spectral transform is performed
using the keyword ``shift-invert`` i.e. the original problem is transformed into
an equivalent problem with eigenvalues given by :math:`\dfrac{1}{\lambda - \sigma}`
instead of :math:`\lambda` where :math:`\sigma` is the value of the spectral shift.
It is therefore much easier to compute eigenvalues close to :math:`\sigma` i.e.
close to :math:`\sigma = 0` in the present case. Eigenvalues are then
transformed back by SLEPc to their original value :math:`\lambda`.


We now ask SLEPc to extract the first 6 eigenvalues by calling its solve function
and extract the corresponding eigenpair (first two arguments of ``get_eigenpair``
correspond to the real and complex part of the eigenvalue, the last two to the
real and complex part of the eigenvector)::

 N_eig = 6   # number of eigenvalues
 print "Computing %i first eigenvalues..." % N_eig
 eigensolver.solve(N_eig)

 # Exact solution computation
 from scipy.optimize import root
 from math import cos, cosh
 falpha = lambda x: cos(x)*cosh(x)+1
 alpha = lambda n: root(falpha, (2*n+1)*pi/2.)['x'][0]

 # Set up file for exporting results
 file_results = XDMFFile("modal_analysis.xdmf")
 file_results.parameters["flush_output"] = True
 file_results.parameters["functions_share_mesh"] = True

 # Extraction
 for i in range(N_eig):
     # Extract eigenpair
     r, c, rx, cx = eigensolver.get_eigenpair(i)

     # 3D eigenfrequency
     freq_3D = sqrt(r)/2/pi

     # Beam eigenfrequency
     if i % 2 == 0: # exact solution should correspond to weak axis bending
         I_bend = H*B**3/12.
     else:          #exact solution should correspond to strong axis bending
         I_bend = B*H**3/12.
     freq_beam = alpha(i/2)**2*sqrt(E*I_bend/(rho*B*H*L**4))/2/pi

     print("Solid FE: {0:8.5f} [Hz]   Beam theory: {1:8.5f} [Hz]".format(freq_3D, freq_beam))

     # Initialize function and assign eigenvector (renormalize by stiffness matrix)
     eigenmode = Function(V,name="Eigenvector "+str(i))
     eigenmode.vector()[:] = rx

The beam analytical solution is obtained using the eigenfrequencies of a clamped
beam in bending given by :math:`\omega_n = \alpha_n^2\sqrt{\dfrac{EI}{\rho S L^4}}`
where :math:`S=BH` is the beam section, :math:`I` the bending inertia and
:math:`\alpha_n` is the solution of the following nonlinear equation:

.. math::
 \cos(\alpha)\cosh(\alpha)+1 = 0

the solution of which can be well approximated by :math:`(2n+1)\pi/2` for :math:`n\geq 3`.
Since the beam possesses two bending axis, each solution to the previous equation is
associated with two frequencies, one with bending along the weak axis (:math:`I=I_{\text{weak}} = HB^3/12`)
and the other along the strong axis (:math:`I=I_{\text{strong}} = BH^3/12`). Since :math:`I_{\text{strong}} = 4I_{\text{weak}}`
for the considered numerical values, the strong axis bending frequency will be twice that corresponsing
to bending along the weak axis. The solution :math:`\alpha_n` are computed using the
``scipy.optimize.root`` function with initial guess given by :math:`(2n+1)\pi/2`.

With ``Nx=400``, we obtain the following comparison between the FE eigenfrequencies
and the beam theory eigenfrequencies :


=====  =============  =================
Mode      Eigenfrequencies
-----  --------------------------------
 #     Solid FE [Hz]   Beam theory [Hz]
=====  =============  =================
  1      2.04991           2.01925
  2      4.04854           4.03850
  3      12.81504         12.65443
  4      25.12717         25.30886
  5      35.74168         35.43277
  6      66.94816         70.86554
=====  =============  =================


