# TC24 same-source ideal and kinetic GPEC result

## Accepted global result

Slurm job `920804` completed both twins with exit `0:0`. Both twins use GPEC
`v1.5.7-581-g2e031ae6`, the same executable pair, equilibrium, coils, kinetic
profile, GPEC deck, PENTRC deck, and numerical controls. The only physics
changes in `dcon.in` are:

```diff
- kin_flag = .false.
- con_flag = .false.
+ kin_flag = .true.
+ con_flag = .true.
```

The ideal case has no native kinetic profile. Its control-surface scalar is
`-1.24265549 N m`, the same small anti-Hermitian residue seen in the archived
ideal calculation. The 408 printed complex control coefficients reproduce the
archived ideal output with relative L2 difference `5.53e-7`. The maximum
component difference is `3.58e-6`.

The kinetic case gives:

- global native GPEC torque: `255.783088 N m`
- native profile endpoint: `255.782098 N m`
- absolute closure error: `0.000990 N m`
- relative closure error: `3.87047e-6`

This passes the global torque gate. It is the self-consistent GPEC and PENTRC
NTV torque for the executed TC24 inputs.

## Local profile status

The 555 native samples span `0.0493157 <= psi_N <= 0.990145`. No smoothing,
interpolation, or fitted factor was used. The differential profile contains
large sign-changing spikes:

- maximum `5.99978e8 N m` at `psi_N = 0.487251`
- minimum `-6.81639e8 N m` at `psi_N = 0.655480`

These locations coincide closely with the ideal `m/n = 4/3` and `5/3`
rational surfaces reported by DCON. Smaller features occur near the outer
`m/3` surfaces. The cumulative profile reaches about `3.09e4 N m` locally
before cancellations leave the accepted `255.782 N m` endpoint.

The plots expose these features on signed logarithmic axes. The global
integral is accepted because it closes against GPEC's separately printed
control scalar. The local spike amplitudes remain diagnostic until a radial
and layer resolution study demonstrates convergence.

## NetCDF profile defect

The independent NetCDF serialization is not a valid oracle at the executed
source commit. `cspline_alloc` allocates torque arrays on indices `0:mstep`,
but `gpout_dw` reshapes the complete arrays into `(mstep,2)` before writing
`T` and `dTdpsi`. The extra index shifts the final real entry into the start of
the imaginary plane and drops data at the end.

The failure is visible without interpretation. The ASCII cumulative endpoint
is `255.78209787 N m`. The NetCDF real plane ends one point earlier at
`254.10304173 N m`, while `255.78209787` appears as the first value of the
next plane.

The accepted CSV and plots use `gpec_dw_profile_n3.out`. This ASCII output
iterates over indices `1:mstep` and closes against the independent global
scalar. The NetCDF `T` and `dTdpsi` variables are rejected for this run. The
defect affects profile serialization, not the completed kinetic solve or the
ASCII and global results.

## Numerical warnings

Kinetic DCON completed normally in `4:10` wall time with 16 OpenMP threads.
It reported displaced-surface search errors, found no kinetic singular
surfaces, and raised IEEE overflow, underflow, and denormal flags. No NaN,
infinity, or maximum-double sentinel occurs in the native outputs. This is
different from the rejected legacy `a430ad2` kinetic attempt, which produced
nonfinite matrix entries and no Euler solution.

Both twins report an unstable free-boundary `n=3` mode. This warning belongs to
the executed case and is preserved in the raw logs.

## Physics interpretation

This is a real finite torque, but its channel is NTV. It confirms that the
kinetic, non-self-adjoint GPEC path has a finite anti-Hermitian response while
the ideal path does not have a physical integrated torque.

It is not the resonant-layer electromagnetic torque from a native MARS
perturbed current crossed with its matching perturbed field. It therefore
does not remove the need for Xingting's native MARS current and exact run
inputs.

The sign in the data and plots is the native GPEC sign. No laboratory CCW sign
is assigned here.

## Evidence

The complete 1.1 GB campaign is immutable at:

```text
/mnt/storage/codex-gpec/tc24-gpec-kinetic-2e031ae6-20260730-r3
```

An `rsync --checksum` comparison against the completed cluster workspace
returned no differences. Exact executable, input, and output hashes are in
[RUN_MANIFEST.json](RUN_MANIFEST.json).
