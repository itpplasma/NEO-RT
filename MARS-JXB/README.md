# Native torque evidence

This directory keeps the compact, reviewable evidence for the TC24 torque
work. Raw GPEC and MARS campaigns remain outside the repository under
`/mnt/storage`.

The current GPEC experiment is a same-source ideal and kinetic comparison.
`stage_gpec_kinetic.sh` creates a unique cluster workspace and
`gpec_kinetic_worker.slurm` builds one executable pair, then runs both input
decks. `analyze_gpec_torque.py` reads native GPEC ASCII outputs and plots them
without interpolation or smoothing.

The kinetic GPEC torque is neoclassical toroidal viscosity. It is not the
resonant-layer electromagnetic torque that MARS reports from its perturbed
current and field. The distinction and the ideal torque null are explained in
[GPEC_IDEAL_KINETIC_TORQUE.md](GPEC_IDEAL_KINETIC_TORQUE.md).

Run a new campaign with:

```sh
./MARS-JXB/stage_gpec_kinetic.sh RUN_TAG
```

After archiving the completed campaign, extract evidence with:

```sh
python MARS-JXB/analyze_gpec_torque.py CAMPAIGN_ROOT OUTPUT_DIR
```

The analysis is strict. It rejects a nonfinite profile, a nonmonotonic native
radial grid, and an ideal case that unexpectedly contains a kinetic torque
profile.
