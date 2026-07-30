# Boozer electromagnetic derivations

## Ampère kernel

`generate_boozer_ampere.f90` asks fortsym to differentiate the real and
quadrature fields associated with

\[
\widehat{\mathbf B}(s)\exp[i(n\phi+m\theta)]
\]

through the chart curl operator.  For coordinates ordered as
`(s, phi, theta)`, Ampère's law gives

\[
\begin{aligned}
\widehat J^s &=
  \frac{i(n\widehat B_\theta-m\widehat B_\phi)}{\mu_0\mathcal J},\\
\widehat J^\phi &=
  \frac{i m\widehat B_s-\partial_s\widehat B_\theta}
       {\mu_0\mathcal J},\\
\widehat J^\theta &=
  \frac{\partial_s\widehat B_\phi-i n\widehat B_s}
       {\mu_0\mathcal J}.
\end{aligned}
\]

Here the magnetic components are covariant, the current components are
contravariant, and \(\mathcal J\) is signed.  The generator emits real and
imaginary parts into `src/generated/boozer_ampere_kernel.f90`.

Regenerate and verify the committed kernel with:

```sh
cd derivations/fortsym
fpm run --profile release
python check_generated.py
```

## Toroidal `j×b` kernel

For coordinates ordered `(s, phi, theta)`, define the Jacobian-weighted
contravariant amplitudes

\[
K^i=\mathcal J J^i,\qquad C^i=\mathcal J B^i.
\]

The coordinate-basis cross product gives the volume-weighted covariant
toroidal force directly:

\[
\mathcal J(\mathbf J\times\mathbf B)_\phi
 =K^\theta C^s-K^sC^\theta .
\]

The vector perturbation contract supplies covariant magnetic components.
Using \(\det g=\mathcal J^2\), the required weighted contravariant field is

\[
C^i=\mathcal J g^{ij}B_j
   =\frac{\operatorname{cof}(g)_{ij}B_j}{\mathcal J}.
\]

`generate_boozer_jxb.f90` constructs the metric cofactors symbolically and
emits both this raise and the complex phase average

\[
\frac12\operatorname{Re}
\left(K^\theta C^{s*}-K^sC^{\theta*}\right)
\]

to `src/generated/boozer_jxb_kernel.f90`.  NEO-RT performs the remaining
angular quadrature. Regenerate this kernel with:

```sh
cd derivations/fortsym
fpm run --profile release generate_boozer_jxb
```

The dependency is pinned to the fortsym revision recorded in the manifest and
generated-file banners. NEO-RT builds the committed leaf kernels and does not
require a computer-algebra package at ordinary build time. Verify both
committed kernels byte-for-byte with:

```sh
python derivations/fortsym/check_generated.py
```
