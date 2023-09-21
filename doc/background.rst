
.. _background:

Background
==========

This page gives a brief review of the discrete Lehmann representation (DLR), and
establishes the definitions and conventions used in ``cppdlr`` (which can vary
from one reference to another). If you are already familiar with the DLR, you
should probably still read :ref:`the section on the relative imaginary time
format below<relativeformat>`. For a more detailed description of the DLR,
please see the references listed on the :ref:`main page<main>`. For examples of
the implementation of these concepts in ``cppdlr``, please see the
:ref:`examples page<examples>`. 

Imaginary time Green's functions and the Lehmann representation
---------------------------------------------------------------

The single-particle imaginary time Green's function is defined in terms of the time-ordered expectation value as

.. math::

   G_{ab}(\tau) = - \langle \mathcal{T} c_a(\tau) c_b^\dagger(0) \rangle, 

where :math:`c^\dagger_b(0)` is the creation operator for a particle in state :math:`b` at time :math:`0` and :math:`c_a(\tau)` is the annihilation operator for a particle in state :math:`a` at time :math:`\tau`. The Green's function is defined on the interval :math:`\tau \in (0, \beta)`, where :math:`\beta` is the inverse temperature, but it can be extended to :math:`\tau \in (-\beta, 0)` using the periodicity and anti-periodicity properties

.. math::
   
   G_{ab}(\tau) = \xi G_{ab}(\beta + \tau),

of bosonic (:math:`\xi = 1`) and fermionic (:math:`\xi = -1`) particles, respectively.

The spectral Lehmann representation of a Green's function is given by

.. math::
   
   G(\tau) = \int_{-\infty}^\infty K(\tau,\omega) \rho(\omega) \, d\omega,

where :math:`\rho(\omega)` is the spectral function corresponding to the Green's
function :math:`G`, and :math:`K(\tau,\omega)` is the analytic continuation
kernel, given by

.. math::

   K(\tau, \omega) = -\frac{e^{-\omega \tau}}{1 + e^{-\beta \omega}}.

Taking the Fourier transform to the imaginary (or Matsubara) frequency domain
gives

.. math::
   
   G(i \nu_n) = \int_{-\infty}^\infty K(i \nu_n,\omega) \rho(\omega) \, d\omega,

with

.. math::

  K(i \nu_n, \omega) = \frac{1}{i\nu_n - \omega}

for fermionic Green's functions, and

.. math::
  K(i \nu_n, \omega) = \frac{\tanh (\beta \omega/2)}{i\nu_n - \omega}

for bosonic Green's functions. Here, the Matsubara frequencies are defined as
:math:`i \nu_n = 2 n \pi i/\beta` and :math:`i \nu_n = (2n+1) \pi i/\beta` for
bosonic and fermionic Green's functions, respectively.

Discrete Lehmann representation
-------------------------------

The discrete Lehmann representation (DLR) is constructed by making a low-rank
approximation of the analytic continuation kernel (using the `interpolative
decomposition <https://epubs.siam.org/doi/10.1137/030602678>`_). Let us define
the dimensionless cutoff paramter :math:`\Lambda \equiv \beta \omega_{\max}`,
where :math:`\omega_{\max}` is defined such that :math:`\rho(\omega) = 0`
outside of :math:`[-\omega_{\max},\omega_{\max}]`. In practice, :math:`\beta` is
typically known, and :math:`\omega_{\max}` can be estimated. :math:`\Lambda` is
a user-specified parameter, and in the typical case that :math:`\omega_{\max}`
is not known exactly, results can be converged with respect to :math:`\Lambda`.
Given :math:`\Lambda` and another user-specified error tolerance parameter
:math:`\epsilon`, the DLR expansion of an imaginary time Green's function
:math:`G(\tau)` is given by

.. math::
  \begin{equation}
    G(\tau) \approx \sum_{l=1}^r K(\tau,\omega_l) \widehat{g}_l, \label{dlrexp} \tag{1}
  \end{equation}

with equality to accuracy :math:`\epsilon`. Here, the :math:`r` *DLR frequencies*
:math:`\omega_l` determining the *DLR basis functions* :math:`K(\tau,\omega_l)`
are carefully chosen (by a pivoted Gram-Schmidt procedure) depending only on
:math:`\Lambda` and :math:`\epsilon`, but *not* on :math:`G` itself. As for the
closely related `intermediate representation
<https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.035147>`_
(implemented, for example, in `sparse-ir
<https://github.com/SpM-lab/sparse-ir>`_), for which the basis functions are
orthogonal but non-explicit, we have :math:`r = \mathcal{O}(\log(\Lambda)
\log(\epsilon^{-1}))`. Thus the DLR enables a highly efficient and high-order
accurate discretization of all imaginary time Green's functions, with a number
of degrees of freedom independent of the specific structure of a given Green's
function (beyond its cutoff parameter :math:`\Lambda`).

Constructing a DLR expansion
----------------------------

In practice, the *DLR coefficients* :math:`\widehat{g}_l` must be determined
from some samples of :math:`G(\tau)`. This can be done by fitting to data, e.g.
via ordinary least squares, or by *interpolation* at the *DLR imaginary time
nodes* :math:`\tau_k`. These are :math:`r` imaginary time points which, given
the DLR frequencies :math:`\omega_l`, are also chosen by a pivoted Gram-Schmidt
procedure. In particular, given the values :math:`G(\tau_k)`, we can solve the
:math:`r \times r` linear system (or interpolation problem)

.. math::
  G(\tau_k) = \sum_{l=1}^r K(\tau_k,\omega_l) \widehat{g}_l

to obtain the DLR coefficients. :math:`G(\tau)` can then be evaluated using its
DLR expansion :math:`\eqref{dlrexp}`.


DLR in the Matsubara frequency domain
-------------------------------------

Fourier transform of :math:`\eqref{dlrexp}` yields

.. math::
  \begin{equation}
    G(i \nu_n) \approx \sum_{l=1}^r K(i \nu_n,\omega_l) \widehat{g}_l, \label{dlrexp_imfreq} \tag{2}
  \end{equation}

so we see that the DLR expansion can be evaluated directly in imaginary time or
imaginary frequency, i.e. the Fourier transform is performed analytically. As in
imaginary time, the DLR coefficients can be obtained by solving the
:math:`r \times r` interpolation problem

.. math::
  G(i \nu_{n_k}) = \sum_{l=1}^r K(i \nu_{n_k},\omega_l) \widehat{g}_l

at the :math:`r` *DLR imaginary frequency nodes* :math:`i \nu_{n_k}`, whereupon
:math:`G(i \nu_n)` can be evaluated using :math:`\eqref{dlrexp_imfreq}` (or
:math:`G(\tau)` can be evaluated using :math:`\eqref{dlrexp}`).


Operations in the DLR basis
---------------------------

Since the DLR basis functions are known analytically, common linear
operations can be straightforwardly performed by representing them in the DLR
basis. These include

- Fourier transform: as explained above, one can switch between imaginary time
  and imaginary frequency representations via the DLR expansion, with no
  additional Fourier transform operation
- Products: in imaginary time or imaginary frequency, by simply multiplying the
  functions on the DLR grid, i.e. :math:`H(\tau_k) = F(\tau_k) G(\tau_k)`,
  whereupon the DLR expansion of the result can be recovered
- Imaginary time convolution: this includes the full convolution
  :math:`H(\tau) = \int_0^\beta F(\tau-\tau') G(\tau') \, d\tau'`, which requires using the
  periodicity/anti-periodicity condition, or the time-ordered convolution
  :math:`H(\tau) = \int_0^\tau F(\tau-\tau') G(\tau') \, d\tau`
- Linear functionals: e.g. inner products with a given function, evaluation at a point, etc...

All such operations take the form of vectors/matrices/tensors acting on :math:`r
\times 1` vectors, which represent the DLR expansion of a Green's function
:math:`G` (either the vector of DLR coefficients of :math:`G`, or the vector of
values of :math:`G` at the DLR nodes). Common operations are implemented in
``cppdlr`` in a user-friendly manner, and the implementation of new operations
should be requested on the `GitHub issues page
<https://github.com/flatironinstitute/cppdlr/issues>`_.


.. _relativeformat:

Imaginary time point format
---------------------------

First, in ``cppdlr`` imaginary time points are scaled from the interval
:math:`[0,\beta]` to the interval :math:`[0,1]`. This is because ``cppdlr``
works with dimensionless variables whenever possible, so in many functions it is
unnecessary to specify the inverse temperature :math:`\beta` explicitly.

Second, ``cppdlr`` stores imaginary time points in a peculiar manner, called the
*relative* time format.
**This is a subtle issue which ``cppdlr`` users should be aware of, in
particular if one wants to supply imaginary time points at which to
evaluate a DLR expansion.** For the TLDR, skip to the **guidelines** below. For an even
more detailed discussion of this issue than the one given here, see Appendix C
of `this paper
<https://www.sciencedirect.com/science/article/pii/S0010465522001771>`_.

The relative time format works as follows. Points :math:`\tau \in [0, 0.5]` are
represented normally. However, instead of representing points :math:`\tau \in
(0.5,1)` directly, we instead store the number :math:`\tau^* = \tau-1`. In other
words, we store the negative distance of :math:`\tau` to 1, rather than tau
itself. Recovering :math:`\tau` in the standard *absolute time format* is
straightforward, and is implemented by the function ``rel2abs``.

The reason for this has to do with maintaining full relative accuracy in
floating point arithmetic. To evaluate the kernel :math:`K(\tau,\omega)`, we
sometimes need to compute the value :math:`1-\tau` for :math:`\tau` very close to 1. If we
work with tau directly, there is a loss of accuracy due to catastrophic
cancellation, which begins to appear in extreme physical regimes and at
very high requested accuracies. If we instead compute :math:`\tau^*` to full relative accuracy and
work with it directly rather than with :math:`\tau`, for example by exploiting
symmetries of :math:`K(\tau,\omega)` to avoid ever evaluating :math:`1-\tau`, we can
maintain full relative accuracy.

This annoyance is the price of maintaining full accuracy in floating point
arithmic. But it is largely ignoreable if the loss of accuracy is not noticeable
in your application, as will be the case for many users.

**Simply follow these guidelines**:

1. Use functions provided by ``cppdlr`` to carry out all imaginary time
   operations whenever possible. This will usually hide this technical
   complication.

2. In a situation in which you want to provide a point :math:`\tau`
   at which to evaluate a DLR, there are two options:

   - (The power user option) Compute :math:`\tau^*`, defined above, to full relative accuracy, and provide this according to
     the instructions in the relevant functions, thereby maintaining full
     relative accuracy in calculations, or
   - (The typical user option) If you don't care about the (usually minor) digit
     loss which comes from ignoring this subtlety, simply convert your point
     :math:`\tau` in the standard, absolute format (a point :math:`\tau \in
     [0,1]`) to the relative format
     :math:`\tau^*` defined above using the ``abs2rel`` function. Since the point will have
     started its life in the absolute format, converting it to relative format
     cannot recover full relative accuracy, but it still needs to be converted
     in order to be compatible with ``cppdlr`` subroutines.

3. If you happen to want to evaluate a Green's function on an
   equispaced grid on :math:`[0,1]` in imaginary time, use the function ``eqpts_rel``
   to generate this grid in the relative format.

Matsubara frequency point format
--------------------------------

We define the Matsubara (or imaginary) frequency points as :math:`i \nu_n = (2 n
+ 1) \pi i/\beta` for fermionic Green's functions, and :math:`i \nu_n = 2 n \pi
i/\beta` for bosonic Green's functions. In ``cppdlr``, Matsubara frequency
points are represented by specifying the integer ``n``, the inverse temperature
:math:`\beta`, and whether the point is a fermionic or bosonic Matsubara
frequency using the ``statistic_t`` specifier.

Symmetrized DLR grids
---------------------

By default, the DLR frequencies :math:`\omega_l` are not chosen to be symmetrized
about :math:`\omega = 0`, nor are the DLR imaginary time nodes :math:`\tau_k`
chosen to be symmetrized about :math:`\tau = \beta/2` or the imaginary frequency
nodes about :math:`i \nu_n = 0`. Indeed, this would represent an additional
constraint in the pivoted Gram-Schmidt procedure used to select the points, and
is not necessary in most applications. However, in some cases, it might be
desirable to have symmetric frequencies and grids, and ``cppdlr`` provides this
functionality via symmetrization flags. Please see the :ref:`list of
other cppdlr capabilities<listofothercapabilities>` section on the :ref:`examples
page<examples>` for a list of `cppdlr` tests which showcase this functionality.

We make a small note about symmetrization for bosonic Green's functions. In this
case, we always select the DLR frequency :math:`\omega = 0`, the DLR imaginary
time node :math:`\tau = \beta/2`, and the DLR imaginary frequency node 
:math:`i \nu_n = 0`. The reason is as follows. A symmetric DLR imaginary frequency grid
containing the point :math:`i \nu_n = 0` must have an odd number of points.
Since it is undesirable to disallow this point, and to make sure all grids have
the same number :math:`r` of points, we force :math:`r` to be odd. This requires
including in each grid, by hand, a single point of high symmetry, and then allowing the
pivoted Gram-Schmidt procedure to select the rest of the points symmetrically.
