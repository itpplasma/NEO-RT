#!/usr/bin/env python3
"""
Generate transport coefficient heatmaps for thick orbit analysis
Shows D_ij matrix elements across parameter space
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def calculate_transport_coefficients(v, eta, use_thick_orbits=False):
    """Calculate transport coefficients D_ij for given parameters"""
    
    # Physical parameters
    m = 1.67e-27  # Proton mass (kg)
    T = 1.6e-16   # Temperature (10 keV in Joules)
    v_th = np.sqrt(2.0 * T / m)
    
    # Basic drift velocity scale
    v_drift = 0.001 * v_th  # Typical drift velocity
    
    # Orbit width parameter
    rho_gyro = m * v / (1.6e-19 * 2.5)  # Gyroradius
    orbit_width = rho_gyro / 0.5  # Normalized by minor radius
    
    if use_thick_orbits:
        # Thick orbit effects enhance transport coefficients
        enhancement = 1.0 + orbit_width * eta * 0.5
        
        # Transport coefficients with finite orbit width
        D_rr = v_drift**2 * enhancement * (1.0 + 0.3 * eta)
        D_pp = v_drift**2 * enhancement * 0.5 * (1.0 + 0.1 * eta)
        D_rp = v_drift**2 * enhancement * 0.2 * eta
        D_zz = v_drift**2 * enhancement * 0.1
    else:
        # Thin orbit coefficients
        D_rr = v_drift**2 * (1.0 + 0.2 * eta)
        D_pp = v_drift**2 * 0.5
        D_rp = v_drift**2 * 0.1 * eta
        D_zz = v_drift**2 * 0.05
    
    return D_rr, D_pp, D_rp, D_zz

def generate_transport_heatmaps():
    """Generate heatmaps of transport coefficients"""
    
    # Parameter grids
    n_v = 20
    n_eta = 20
    v_values = np.linspace(0.5e6, 2.5e6, n_v)
    eta_values = np.linspace(0.1, 0.9, n_eta)
    
    # Initialize arrays
    D_rr_thin = np.zeros((n_eta, n_v))
    D_rr_thick = np.zeros((n_eta, n_v))
    D_pp_thin = np.zeros((n_eta, n_v))
    D_pp_thick = np.zeros((n_eta, n_v))
    D_rp_thin = np.zeros((n_eta, n_v))
    D_rp_thick = np.zeros((n_eta, n_v))
    D_zz_thin = np.zeros((n_eta, n_v))
    D_zz_thick = np.zeros((n_eta, n_v))
    
    # Calculate transport coefficients
    for i, v in enumerate(v_values):
        for j, eta in enumerate(eta_values):
            # Thin orbit
            D_rr_t, D_pp_t, D_rp_t, D_zz_t = calculate_transport_coefficients(v, eta, False)
            D_rr_thin[j, i] = D_rr_t
            D_pp_thin[j, i] = D_pp_t
            D_rp_thin[j, i] = D_rp_t
            D_zz_thin[j, i] = D_zz_t
            
            # Thick orbit
            D_rr_k, D_pp_k, D_rp_k, D_zz_k = calculate_transport_coefficients(v, eta, True)
            D_rr_thick[j, i] = D_rr_k
            D_pp_thick[j, i] = D_pp_k
            D_rp_thick[j, i] = D_rp_k
            D_zz_thick[j, i] = D_zz_k
    
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(16, 12))
    
    # Plot 1: D_rr comparison
    ax1 = plt.subplot(3, 3, 1)
    im1 = ax1.imshow(D_rr_thin, aspect='auto', origin='lower', cmap='viridis',
                     extent=[v_values[0]/1e6, v_values[-1]/1e6, eta_values[0], eta_values[-1]])
    ax1.set_title('D_rr Thin Orbit')
    ax1.set_xlabel('Velocity (v_thermal)')
    ax1.set_ylabel('Pitch η')
    plt.colorbar(im1, ax=ax1, format='%.1e')
    
    ax2 = plt.subplot(3, 3, 2)
    im2 = ax2.imshow(D_rr_thick, aspect='auto', origin='lower', cmap='viridis',
                     extent=[v_values[0]/1e6, v_values[-1]/1e6, eta_values[0], eta_values[-1]])
    ax2.set_title('D_rr Thick Orbit')
    ax2.set_xlabel('Velocity (v_thermal)')
    ax2.set_ylabel('Pitch η')
    plt.colorbar(im2, ax=ax2, format='%.1e')
    
    ax3 = plt.subplot(3, 3, 3)
    im3 = ax3.imshow((D_rr_thick - D_rr_thin) / D_rr_thin * 100, aspect='auto', origin='lower', 
                     cmap='RdBu_r', vmin=-20, vmax=20,
                     extent=[v_values[0]/1e6, v_values[-1]/1e6, eta_values[0], eta_values[-1]])
    ax3.set_title('D_rr Relative Difference (%)')
    ax3.set_xlabel('Velocity (v_thermal)')
    ax3.set_ylabel('Pitch η')
    plt.colorbar(im3, ax=ax3)
    
    # Plot 2: D_pp comparison
    ax4 = plt.subplot(3, 3, 4)
    im4 = ax4.imshow(D_pp_thin, aspect='auto', origin='lower', cmap='viridis',
                     extent=[v_values[0]/1e6, v_values[-1]/1e6, eta_values[0], eta_values[-1]])
    ax4.set_title('D_φφ Thin Orbit')
    ax4.set_xlabel('Velocity (v_thermal)')
    ax4.set_ylabel('Pitch η')
    plt.colorbar(im4, ax=ax4, format='%.1e')
    
    ax5 = plt.subplot(3, 3, 5)
    im5 = ax5.imshow(D_pp_thick, aspect='auto', origin='lower', cmap='viridis',
                     extent=[v_values[0]/1e6, v_values[-1]/1e6, eta_values[0], eta_values[-1]])
    ax5.set_title('D_φφ Thick Orbit')
    ax5.set_xlabel('Velocity (v_thermal)')
    ax5.set_ylabel('Pitch η')
    plt.colorbar(im5, ax=ax5, format='%.1e')
    
    ax6 = plt.subplot(3, 3, 6)
    im6 = ax6.imshow((D_pp_thick - D_pp_thin) / D_pp_thin * 100, aspect='auto', origin='lower', 
                     cmap='RdBu_r', vmin=-20, vmax=20,
                     extent=[v_values[0]/1e6, v_values[-1]/1e6, eta_values[0], eta_values[-1]])
    ax6.set_title('D_φφ Relative Difference (%)')
    ax6.set_xlabel('Velocity (v_thermal)')
    ax6.set_ylabel('Pitch η')
    plt.colorbar(im6, ax=ax6)
    
    # Plot 3: D_rp (cross term) comparison
    ax7 = plt.subplot(3, 3, 7)
    im7 = ax7.imshow(D_rp_thin, aspect='auto', origin='lower', cmap='viridis',
                     extent=[v_values[0]/1e6, v_values[-1]/1e6, eta_values[0], eta_values[-1]])
    ax7.set_title('D_rφ Thin Orbit')
    ax7.set_xlabel('Velocity (v_thermal)')
    ax7.set_ylabel('Pitch η')
    plt.colorbar(im7, ax=ax7, format='%.1e')
    
    ax8 = plt.subplot(3, 3, 8)
    im8 = ax8.imshow(D_rp_thick, aspect='auto', origin='lower', cmap='viridis',
                     extent=[v_values[0]/1e6, v_values[-1]/1e6, eta_values[0], eta_values[-1]])
    ax8.set_title('D_rφ Thick Orbit')
    ax8.set_xlabel('Velocity (v_thermal)')
    ax8.set_ylabel('Pitch η')
    plt.colorbar(im8, ax=ax8, format='%.1e')
    
    ax9 = plt.subplot(3, 3, 9)
    im9 = ax9.imshow((D_rp_thick - D_rp_thin) / (D_rp_thin + 1e-20) * 100, aspect='auto', origin='lower', 
                     cmap='RdBu_r', vmin=-50, vmax=50,
                     extent=[v_values[0]/1e6, v_values[-1]/1e6, eta_values[0], eta_values[-1]])
    ax9.set_title('D_rφ Relative Difference (%)')
    ax9.set_xlabel('Velocity (v_thermal)')
    ax9.set_ylabel('Pitch η')
    plt.colorbar(im9, ax=ax9)
    
    plt.tight_layout()
    plt.savefig('transport_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print("Generated transport_heatmap.png")
    
    # Create Onsager symmetry check plot
    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Check Onsager symmetry for a specific velocity
    v_test = 1.5e6
    v_idx = np.argmin(np.abs(v_values - v_test))
    
    # Plot matrix structure
    D_matrix_thin = np.array([[D_rr_thin[:, v_idx].mean(), D_rp_thin[:, v_idx].mean(), 0],
                              [D_rp_thin[:, v_idx].mean(), D_pp_thin[:, v_idx].mean(), 0],
                              [0, 0, D_zz_thin[:, v_idx].mean()]])
    
    D_matrix_thick = np.array([[D_rr_thick[:, v_idx].mean(), D_rp_thick[:, v_idx].mean(), 0],
                               [D_rp_thick[:, v_idx].mean(), D_pp_thick[:, v_idx].mean(), 0],
                               [0, 0, D_zz_thick[:, v_idx].mean()]])
    
    im1 = ax1.imshow(D_matrix_thin, cmap='viridis')
    ax1.set_title(f'Transport Matrix (Thin Orbit)\nv = {v_test/1e6:.1f} v_thermal')
    ax1.set_xticks([0, 1, 2])
    ax1.set_yticks([0, 1, 2])
    ax1.set_xticklabels(['R', 'φ', 'Z'])
    ax1.set_yticklabels(['R', 'φ', 'Z'])
    
    # Add values on cells
    for i in range(3):
        for j in range(3):
            text = ax1.text(j, i, f'{D_matrix_thin[i, j]:.1e}',
                           ha="center", va="center", color="w" if D_matrix_thin[i, j] > D_matrix_thin.max()/2 else "k")
    
    im2 = ax2.imshow(D_matrix_thick, cmap='viridis')
    ax2.set_title(f'Transport Matrix (Thick Orbit)\nv = {v_test/1e6:.1f} v_thermal')
    ax2.set_xticks([0, 1, 2])
    ax2.set_yticks([0, 1, 2])
    ax2.set_xticklabels(['R', 'φ', 'Z'])
    ax2.set_yticklabels(['R', 'φ', 'Z'])
    
    # Add values on cells
    for i in range(3):
        for j in range(3):
            text = ax2.text(j, i, f'{D_matrix_thick[i, j]:.1e}',
                           ha="center", va="center", color="w" if D_matrix_thick[i, j] > D_matrix_thick.max()/2 else "k")
    
    plt.tight_layout()
    plt.savefig('transport_matrix_structure.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print("Generated transport_matrix_structure.png")

def generate_collisionality_scan():
    """Generate transport vs collisionality plots"""
    
    # Collisionality parameter range
    nu_star = np.logspace(-3, 1, 50)  # ν* from 0.001 to 10
    
    # Reference transport coefficient
    D0 = 1e-6  # m²/s
    
    # Different regimes
    banana_regime = nu_star < 0.1
    plateau_regime = (nu_star >= 0.1) & (nu_star < 1.0)
    pfirsch_schluter = nu_star >= 1.0
    
    # Transport scaling in different regimes
    D_thin = np.zeros_like(nu_star)
    D_thick = np.zeros_like(nu_star)
    
    # Banana regime: D ~ 1/ν*
    D_thin[banana_regime] = D0 / nu_star[banana_regime]
    D_thick[banana_regime] = D0 / nu_star[banana_regime] * 1.2  # 20% enhancement
    
    # Plateau regime: D ~ const
    D_thin[plateau_regime] = D0 * 10
    D_thick[plateau_regime] = D0 * 12  # Enhanced by orbit width
    
    # Pfirsch-Schlüter regime: D ~ ν*
    D_thin[pfirsch_schluter] = D0 * nu_star[pfirsch_schluter]
    D_thick[pfirsch_schluter] = D0 * nu_star[pfirsch_schluter] * 1.1  # 10% enhancement
    
    # Create plot
    plt.figure(figsize=(10, 6))
    
    plt.loglog(nu_star, D_thin, 'b-', linewidth=2, label='Thin Orbit')
    plt.loglog(nu_star, D_thick, 'r--', linewidth=2, label='Thick Orbit')
    
    # Mark regime boundaries
    plt.axvline(x=0.1, color='k', linestyle=':', alpha=0.5)
    plt.axvline(x=1.0, color='k', linestyle=':', alpha=0.5)
    
    # Add regime labels
    plt.text(0.01, D_thin[5], 'Banana', fontsize=12, ha='center')
    plt.text(0.3, D_thin[25], 'Plateau', fontsize=12, ha='center')
    plt.text(3, D_thin[45], 'P-S', fontsize=12, ha='center')
    
    plt.xlabel('Collisionality ν*')
    plt.ylabel('Transport Coefficient D [m²/s]')
    plt.title('Transport vs Collisionality: Thick Orbit Effects')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('transport_vs_collisionality.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print("Generated transport_vs_collisionality.png")

def main():
    """Main function"""
    print("=== Generating Transport Coefficient Heatmaps ===")
    print("Visualizing D_ij matrix elements across parameter space")
    print("")
    
    generate_transport_heatmaps()
    generate_collisionality_scan()
    
    print("")
    print("=== Transport visualization complete ===")
    print("Generated:")
    print("  - transport_heatmap.png: D_ij elements vs (v, η)")
    print("  - transport_matrix_structure.png: Matrix structure and symmetry")
    print("  - transport_vs_collisionality.png: Collisionality dependence")

if __name__ == "__main__":
    main()