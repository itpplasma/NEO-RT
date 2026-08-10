#!/usr/bin/env python3
"""
Generate POTATO profile_poly.in from NEO-2 HDF5 output.

Reads neo2_out.h5 and computes the electrostatic potential phi_e from the
ExB rotation (MtOvR):

    Er(s) = sign_theta * psi_pr * vth(s) * MtOvR_ion(s) / C_CGS

where psi_pr must be the POLOIDAL flux difference in POTATO's own gauge,
psi_pr = (PsiedgeVs - PsiaxisVs) * 1e8 Mx from the EQDSK header, because
POTATO's phielec_of_psi divides dPhi/ds_pol by exactly that quantity.
Using NEO-2's psi_pr_hat*Bref here instead (that is the TOROIDAL flux
derivative, -3.585e7 Mx) inflates Omega_ExB by x1.945 at every radius --
the bug that made the AUG #30835 POTATO/NEO-RT torque benchmark disagree
by 1.6x until 2026-07-16.  NB neort_to_potato.py still has this bug, and
additionally fits in s_tor without converting to s_pol -- prefer this
script.

Usage:
    python import_neo2_profiles.py [--neo2 neo2_out.h5] [--sign-theta -1]
"""

import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt
from libneo import FluxConverter, eqdsk_base

C_CGS = 2.9979e10      # cm/s
POLY_ORDER = 9
SPECIES_ION = 1        # 0-indexed: 0=electrons, 1=ions


def read_neo2(path):
    with h5py.File(path, 'r') as f:
        s_tor = f['boozer_s'][:]
        aiota = f['aiota'][:]
        n_spec = f['n_spec'][:]          # (n_s, 2), 1/cm^3
        T_spec = f['T_spec'][:]          # (n_s, 2), erg
        m_spec = f['m_spec'][:]          # (n_s, 2), g
        MtOvR = f['MtOvR'][:]           # (n_s, 2), 1/cm
        Bref = f['Bref'][:]              # (n_s,), G
        psi_pr_hat = f['psi_pr_hat'][:]  # (n_s,), cm^2  (= psi_pr / Bref)
        R0 = f['R0'][0]                  # cm, major radius (constant flux function)
    return s_tor, aiota, n_spec, T_spec, m_spec, MtOvR, Bref, psi_pr_hat, R0


def stor_to_spol(s_tor, eqdsk_path):
    """Convert normalized toroidal flux to normalized poloidal flux using the
    libneo FluxConverter built from the EFIT q-profile, so the profile radial
    coordinate is consistent with POTATO's own field equilibrium.

    The previous hand-rolled integral used ds_pol ~ q ds_tor, which is inverted
    (the correct relation is ds_pol ~ aiota ds_tor = ds_tor / q); that bug
    radially misplaced every profile (e.g. the Er=0 surface from rho_pol~0.88
    to ~0.72).
    """
    eq = eqdsk_base.read_eqdsk(eqdsk_path)
    converter = FluxConverter(eq["qprof"], eq["s_pol"])
    return np.asarray(converter.stor2spol(s_tor))


def compute_er(T_ion, m_ion, MtOvR_ion, psi_pr, sign_theta):
    vth = np.sqrt(2.0 * T_ion / m_ion)
    return sign_theta * psi_pr * vth * MtOvR_ion / C_CGS


def integrate_er(s, er):
    idx = np.argsort(s)
    s_d, er_d = s[idx], er[idx]
    phi = np.zeros(len(s))
    for i in range(1, len(s)):
        phi[i] = phi[i - 1] - 0.5 * (er_d[i - 1] + er_d[i]) * (s_d[i] - s_d[i - 1])
    result = np.empty(len(s))
    result[idx] = phi
    return result


def polyfit_potato(s, y, order):
    return np.polyfit(s, y, order)


def write_profile_poly(path, poly_n, poly_Te, poly_Ti, poly_phi, order):
    def fmt(coeffs):
        return '  ' + '  '.join(f'{c: .15e}' for c in coeffs)
    with open(path, 'w') as f:
        f.write(' # Generated from NEO-2 HDF5 by import_neo2_profiles.py\n')
        f.write(f' # Polynomial order: {order}  phi_e from MtOvR (ExB rotation)\n')
        f.write(fmt(poly_n) + '\n')
        f.write(fmt(poly_Te) + '\n')
        f.write(fmt(poly_Ti) + '\n')
        f.write(fmt(poly_phi) + '\n')


def save_plots(spol, n_ion, T_ion_eV, Te_eV, er, phi_e,
               poly_n, poly_Ti, poly_phi, order):
    s_fine = np.linspace(0.0, 1.0, 300)

    fig, ax = plt.subplots()
    ax.plot(spol, n_ion, 'ro', ms=3, label='NEO-2 data')
    ax.plot(s_fine, np.polyval(poly_n, s_fine), 'b-', label='polynomial fit')
    ax.set_xlabel('Normalised poloidal flux')
    ax.set_ylabel('Ion density [cm⁻³]')
    ax.set_title('Ion density profile')
    ax.legend()
    fig.savefig('density_fit.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(spol, T_ion_eV, 'ro', ms=3, label='NEO-2 data')
    ax.plot(s_fine, np.polyval(poly_Ti, s_fine), 'b-', label='polynomial fit')
    ax.set_xlabel('Normalised poloidal flux')
    ax.set_ylabel('Ion temperature [eV]')
    ax.set_title('Ion temperature profile')
    ax.legend()
    fig.savefig('temperature_fit.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    er_poly = np.polyfit(spol, er, order + 1)
    fig, ax = plt.subplots()
    ax.plot(spol, er, 'ro', ms=3, label='from MtOvR')
    ax.plot(s_fine, np.polyval(er_poly, s_fine), 'b-', label='polynomial fit')
    ax.set_xlabel('Normalised poloidal flux')
    ax.set_ylabel('Er [statV per unit s]')
    ax.set_title('Radial electric field (from ExB rotation)')
    ax.legend()
    fig.savefig('Er_fit.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(spol, phi_e, 'ro', ms=3, label='integrated Er')
    ax.plot(s_fine, np.polyval(poly_phi, s_fine), 'b-', label='polynomial fit')
    ax.set_xlabel('Normalised poloidal flux')
    ax.set_ylabel('φₑ [statV]')
    ax.set_title('Electrostatic potential')
    ax.legend()
    fig.savefig('phi_e_fit.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    print('Saved density_fit.png, temperature_fit.png, Er_fit.png, phi_e_fit.png')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--neo2', default='neo2_out.h5')
    ap.add_argument('--eqdsk', default='efit.eqdsk',
                    help='EFIT g-file; s_tor<->s_pol via libneo FluxConverter '
                         '(matches POTATO field equilibrium)')
    ap.add_argument('--sign-theta', type=float, default=-1.0,
                    help='Sign of poloidal field (default: -1, AUG standard)')
    ap.add_argument('--output', default='profile_poly.in')
    ap.add_argument('--shift_mt', type=float, default=0.0,
                    help='Constant M_t shift (unscaled) to remove Er=0 crossing')
    ap.add_argument('--plot', action='store_true')
    args = ap.parse_args()

    s_tor, aiota, n_spec, T_spec, m_spec, MtOvR, Bref, psi_pr_hat, R0 = read_neo2(args.neo2)
    spol = stor_to_spol(s_tor, args.eqdsk)

    n_ion = n_spec[:, SPECIES_ION]
    T_ion = T_spec[:, SPECIES_ION]   # erg
    T_ion_eV = T_ion / 1.6e-12
    Te = T_spec[:, 0]                 # erg
    Te_eV = Te / 1.6e-12
    m_ion = m_spec[:, SPECIES_ION]
    MtOvR_ion = MtOvR[:, SPECIES_ION] + args.shift_mt / R0  # apply constant M_t shift

    # psi_pr must be the SAME poloidal-flux normalization POTATO divides by in
    # phielec_of_psi (field_eq psi_sep - psi_axis, EQDSK gauge, Mx): using
    # NEO-2's psi_pr_hat[0]*Bref[0] here inflated Er by x1.945 uniformly.
    eq_hdr = eqdsk_base.read_eqdsk(args.eqdsk)
    psi_pr = (eq_hdr['PsiedgeVs'] - eq_hdr['PsiaxisVs']) * 1.0e8  # Mx
    er = compute_er(T_ion, m_ion, MtOvR_ion, psi_pr, args.sign_theta)
    phi_e = integrate_er(spol, er)

    poly_n = polyfit_potato(spol, n_ion, POLY_ORDER)
    poly_Ti = polyfit_potato(spol, T_ion_eV, POLY_ORDER)
    poly_Te = polyfit_potato(spol, Te_eV, POLY_ORDER)
    poly_phi = polyfit_potato(spol, phi_e, POLY_ORDER)

    write_profile_poly(args.output, poly_n, poly_Te, poly_Ti, poly_phi, POLY_ORDER)
    Mt = MtOvR_ion * R0
    print(f'Written {args.output}')
    print(f'M_t shift = {args.shift_mt:+.4e}, R0 = {R0:.2f} cm')
    print(f'M_t range: [{Mt.min():+.4e}, {Mt.max():+.4e}] (zero crossing: {(Mt.min()<0)*(Mt.max()>0)})')
    print(f'phi_e range: {phi_e.min():.4f} to {phi_e.max():.4f} statV')
    print(f'Er range: {er.min():.4f} to {er.max():.4f} statV/unit_s')
    print(f'psi_pr = {psi_pr:.4e} Mx, sign_theta = {args.sign_theta:+.0f}')

    if args.plot:
        save_plots(spol, n_ion, T_ion_eV, Te_eV, er, phi_e,
                   poly_n, poly_Ti, poly_phi, POLY_ORDER)


if __name__ == '__main__':
    main()
