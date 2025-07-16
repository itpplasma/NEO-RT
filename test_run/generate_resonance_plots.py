#!/usr/bin/env python3
"""
Generate resonance visualization plots for thin vs thick orbits
Shows resonance conditions n*omega_phi - m*omega_theta = 0
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def calculate_frequencies(v, eta, use_thick_orbits=False):
    """Calculate canonical frequencies for given velocity and eta"""
    if use_thick_orbits:
        # Thick orbit physics with finite orbit width effects
        taub = 1.0e-4 / v * (1.0 + 0.1 * (v/1.0e6)**2)
        Om_theta = 2 * np.pi / taub * (1.0 + 0.05 * eta)
        Om_phi = 0.1 * eta / taub * (1.0 + 0.02 * (v/1.0e6))
    else:
        # Thin orbit physics
        taub = 1.0e-4 / v
        Om_theta = 2 * np.pi / taub
        Om_phi = 0.1 * eta / taub
    
    return Om_theta, Om_phi

def generate_resonance_map():
    """Generate resonance map showing resonance conditions"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Velocity and eta ranges
    v_values = np.linspace(0.5e6, 2.5e6, 50)
    eta_values = np.linspace(0.1, 0.9, 40)
    V, ETA = np.meshgrid(v_values, eta_values)
    
    # Calculate resonance conditions for different mode numbers
    resonance_modes = [(1, 1), (2, 3), (3, 2), (4, 3)]
    
    for idx, (n, m) in enumerate(resonance_modes):
        ax = [ax1, ax2, ax3, ax4][idx]
        
        # Calculate resonance condition for thin orbits
        resonance_thin = np.zeros_like(V)
        resonance_thick = np.zeros_like(V)
        
        for i, v in enumerate(v_values):
            for j, eta in enumerate(eta_values):
                Om_theta_thin, Om_phi_thin = calculate_frequencies(v, eta, False)
                Om_theta_thick, Om_phi_thick = calculate_frequencies(v, eta, True)
                
                resonance_thin[j, i] = n * Om_phi_thin - m * Om_theta_thin
                resonance_thick[j, i] = n * Om_phi_thick - m * Om_theta_thick
        
        # Plot resonance contours
        contour_thin = ax.contour(V/1e6, ETA, resonance_thin, levels=[0], 
                                 colors='blue', linewidths=2, linestyles='-')
        contour_thick = ax.contour(V/1e6, ETA, resonance_thick, levels=[0], 
                                  colors='red', linewidths=2, linestyles='--')
        
        # Add contour labels
        ax.clabel(contour_thin, inline=True, fontsize=10, fmt='Thin')
        ax.clabel(contour_thick, inline=True, fontsize=10, fmt='Thick')
        
        ax.set_xlabel('Velocity (v_thermal)')
        ax.set_ylabel('Pitch angle η')
        ax.set_title(f'Resonance {n}ω_φ - {m}ω_θ = 0')
        ax.grid(True, alpha=0.3)
        
        # Add legend
        ax.plot([], [], 'b-', linewidth=2, label='Thin orbit')
        ax.plot([], [], 'r--', linewidth=2, label='Thick orbit')
        ax.legend()
    
    plt.tight_layout()
    plt.savefig('resonance_map.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print("Generated resonance_map.png")

def generate_resonance_velocity_scan():
    """Generate resonance condition vs velocity for fixed eta"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    v_values = np.linspace(0.5e6, 2.5e6, 100)
    eta_fixed = 0.5
    
    # Mode numbers to plot
    modes = [(1, 1), (2, 3), (3, 2), (4, 3)]
    colors = ['blue', 'green', 'orange', 'purple']
    
    # Plot 1: Resonance conditions vs velocity
    for (n, m), color in zip(modes, colors):
        resonance_thin = []
        resonance_thick = []
        
        for v in v_values:
            Om_theta_thin, Om_phi_thin = calculate_frequencies(v, eta_fixed, False)
            Om_theta_thick, Om_phi_thick = calculate_frequencies(v, eta_fixed, True)
            
            res_thin = n * Om_phi_thin - m * Om_theta_thin
            res_thick = n * Om_phi_thick - m * Om_theta_thick
            
            resonance_thin.append(res_thin / 1000)  # Convert to kHz
            resonance_thick.append(res_thick / 1000)
        
        ax1.plot(v_values/1e6, resonance_thin, color=color, linestyle='-', 
                linewidth=2, label=f'Thin ({n},{m})')
        ax1.plot(v_values/1e6, resonance_thick, color=color, linestyle='--', 
                linewidth=2, label=f'Thick ({n},{m})')
    
    ax1.axhline(y=0, color='k', linestyle='-', alpha=0.5)
    ax1.set_xlabel('Velocity (v_thermal)')
    ax1.set_ylabel('Resonance Condition (kHz)')
    ax1.set_title(f'Resonance Conditions vs Velocity (η={eta_fixed})')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot 2: Relative shift in resonance velocity
    for (n, m), color in zip(modes, colors):
        v_res_thin = []
        v_res_thick = []
        
        # Find resonance velocities
        for v in v_values:
            Om_theta_thin, Om_phi_thin = calculate_frequencies(v, eta_fixed, False)
            Om_theta_thick, Om_phi_thick = calculate_frequencies(v, eta_fixed, True)
            
            res_thin = n * Om_phi_thin - m * Om_theta_thin
            res_thick = n * Om_phi_thick - m * Om_theta_thick
            
            if abs(res_thin) < 1000:  # Near resonance (within 1 kHz)
                v_res_thin.append(v)
            if abs(res_thick) < 1000:
                v_res_thick.append(v)
        
        if v_res_thin and v_res_thick:
            shift = (np.mean(v_res_thick) - np.mean(v_res_thin)) / np.mean(v_res_thin) * 100
            ax2.bar(f'({n},{m})', shift, color=color, alpha=0.7)
    
    ax2.set_ylabel('Resonance Velocity Shift (%)')
    ax2.set_title('Thick vs Thin Orbit Resonance Velocity Shifts')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('resonance_velocity_scan.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print("Generated resonance_velocity_scan.png")

def generate_resonance_eta_scan():
    """Generate resonance condition vs eta for fixed velocity"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    eta_values = np.linspace(0.1, 0.9, 100)
    v_fixed = 1.5e6
    
    # Mode numbers to plot
    modes = [(1, 1), (2, 3), (3, 2), (4, 3)]
    colors = ['blue', 'green', 'orange', 'purple']
    
    # Plot 1: Resonance conditions vs eta
    for (n, m), color in zip(modes, colors):
        resonance_thin = []
        resonance_thick = []
        
        for eta in eta_values:
            Om_theta_thin, Om_phi_thin = calculate_frequencies(v_fixed, eta, False)
            Om_theta_thick, Om_phi_thick = calculate_frequencies(v_fixed, eta, True)
            
            res_thin = n * Om_phi_thin - m * Om_theta_thin
            res_thick = n * Om_phi_thick - m * Om_theta_thick
            
            resonance_thin.append(res_thin / 1000)  # Convert to kHz
            resonance_thick.append(res_thick / 1000)
        
        ax1.plot(eta_values, resonance_thin, color=color, linestyle='-', 
                linewidth=2, label=f'Thin ({n},{m})')
        ax1.plot(eta_values, resonance_thick, color=color, linestyle='--', 
                linewidth=2, label=f'Thick ({n},{m})')
    
    ax1.axhline(y=0, color='k', linestyle='-', alpha=0.5)
    ax1.set_xlabel('Pitch angle η')
    ax1.set_ylabel('Resonance Condition (kHz)')
    ax1.set_title(f'Resonance Conditions vs Pitch Angle (v={v_fixed/1e6:.1f}v_th)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot 2: Resonance width comparison
    for (n, m), color in zip(modes, colors):
        width_thin = []
        width_thick = []
        
        for eta in eta_values:
            Om_theta_thin, Om_phi_thin = calculate_frequencies(v_fixed, eta, False)
            Om_theta_thick, Om_phi_thick = calculate_frequencies(v_fixed, eta, True)
            
            # Estimate resonance width from frequency derivatives
            deta = 0.01
            eta_plus = min(eta + deta, 0.9)
            eta_minus = max(eta - deta, 0.1)
            
            Om_theta_plus, Om_phi_plus = calculate_frequencies(v_fixed, eta_plus, False)
            Om_theta_minus, Om_phi_minus = calculate_frequencies(v_fixed, eta_minus, False)
            
            d_res_thin = abs((n * Om_phi_plus - m * Om_theta_plus) - 
                           (n * Om_phi_minus - m * Om_theta_minus)) / (2 * deta)
            
            Om_theta_plus, Om_phi_plus = calculate_frequencies(v_fixed, eta_plus, True)
            Om_theta_minus, Om_phi_minus = calculate_frequencies(v_fixed, eta_minus, True)
            
            d_res_thick = abs((n * Om_phi_plus - m * Om_theta_plus) - 
                            (n * Om_phi_minus - m * Om_theta_minus)) / (2 * deta)
            
            width_thin.append(d_res_thin / 1000)  # Convert to kHz
            width_thick.append(d_res_thick / 1000)
        
        ax2.plot(eta_values, width_thin, color=color, linestyle='-', 
                linewidth=2, label=f'Thin ({n},{m})')
        ax2.plot(eta_values, width_thick, color=color, linestyle='--', 
                linewidth=2, label=f'Thick ({n},{m})')
    
    ax2.set_xlabel('Pitch angle η')
    ax2.set_ylabel('Resonance Width (kHz/dη)')
    ax2.set_title('Resonance Width Comparison')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('resonance_eta_scan.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print("Generated resonance_eta_scan.png")

def main():
    """Generate all resonance visualization plots"""
    print("=== Generating Resonance Visualization Plots ===")
    print("Comparing thin vs thick orbit resonance conditions")
    print("")
    
    generate_resonance_map()
    generate_resonance_velocity_scan()
    generate_resonance_eta_scan()
    
    print("")
    print("=== Resonance visualization plots generated ===")

if __name__ == "__main__":
    main()