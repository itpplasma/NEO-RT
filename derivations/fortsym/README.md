# Boozer Ampère derivation

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

The dependency is pinned to the fortsym revision recorded in both the manifest
and generated-file banner.  NEO-RT builds the committed leaf kernel and does
not require a computer-algebra package at ordinary build time.
