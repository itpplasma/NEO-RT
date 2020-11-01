Code design
===========

Bounce averages are required to compute canonical frequencies
and Fourier harmonics. This is implemented in :f:subr:`~orbit/bounce_average`
and relies on a given starting position :code:`zstart`.
A suitable starting position with given invariants :math:`\balpha`
is provided by :code:`TODO`.

1. Read fields and profiles
2. Precompute and interpolate canonical frequencies :math:`\Omega^k(\balpha)`