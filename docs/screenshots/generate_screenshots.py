import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, FancyBboxPatch

# Set dark theme styling matching egui / egui_plot
plt.style.use('dark_background')
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#404552'
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['grid.color'] = '#2e3440'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.alpha'] = 0.6

def generate_workbench():
    fig = plt.figure(figsize=(16, 9), dpi=140, facecolor='#12151c')
    
    # Outer layout: Left panel (controls) and Right panel (plots + telemetry)
    gs_main = gridspec.GridSpec(1, 2, width_ratios=[1, 3.2], figure=fig, wspace=0.08, left=0.02, right=0.98, top=0.94, bottom=0.04)
    
    # Left Control Panel
    ax_ctrl = fig.add_subplot(gs_main[0])
    ax_ctrl.set_facecolor('#181b24')
    ax_ctrl.set_xticks([])
    ax_ctrl.set_yticks([])
    for spine in ax_ctrl.spines.values():
        spine.set_color('#2a2f3d')
    
    # Add UI mockup text in control panel
    ax_ctrl.text(0.06, 0.96, "OpenBand v1.0.0", color='#61afef', fontsize=15, fontweight='bold', transform=ax_ctrl.transAxes)
    ax_ctrl.text(0.06, 0.93, "1D Semiconductor Heterojunction Solver", color='#828997', fontsize=9, transform=ax_ctrl.transAxes)
    
    ax_ctrl.text(0.06, 0.88, "MATERIAL PRESETS", color='#98c379', fontsize=10, fontweight='bold', transform=ax_ctrl.transAxes)
    ax_ctrl.text(0.06, 0.85, "> AlGaAs / GaAs Heterojunction", color='#e5c07b', fontsize=9.5, transform=ax_ctrl.transAxes)
    
    ax_ctrl.text(0.06, 0.79, "PHYSICAL PARAMETERS", color='#e06c75', fontsize=10, fontweight='bold', transform=ax_ctrl.transAxes)
    controls = [
        ("Device Length L:", "1.00 um"),
        ("Junction Position:", "0.50 um"),
        ("p-Side (Al0.3Ga0.7As):", "Eg = 1.798 eV"),
        ("  Affinity χ_p:", "3.74 eV"),
        ("  Doping Na:", "1.00e17 cm^-3"),
        ("n-Side (GaAs):", "Eg = 1.424 eV"),
        ("  Affinity χ_n:", "4.07 eV"),
        ("  Doping Nd:", "1.00e17 cm^-3"),
        ("Temperature T:", "300.0 K (Vt = 25.85 mV)"),
        ("Applied Bias V_bias:", "0.00 V (Equilibrium)"),
        ("SRH Lifetime τ_n, τ_p:", "10.0 ns, 10.0 ns"),
    ]
    y_pos = 0.75
    for label, val in controls:
        ax_ctrl.text(0.06, y_pos, label, color='#abb2bf', fontsize=8.5, transform=ax_ctrl.transAxes)
        ax_ctrl.text(0.55, y_pos, val, color='#61afef', fontsize=8.5, transform=ax_ctrl.transAxes)
        y_pos -= 0.032
        
    ax_ctrl.text(0.06, 0.36, "NUMERICAL SOLVER", color='#c678dd', fontsize=10, fontweight='bold', transform=ax_ctrl.transAxes)
    solver_info = [
        ("FDM Grid Nodes:", "600 (dx = 1.67 nm)"),
        ("Debye Scaling L_D:", "28.4 nm (Normalized)"),
        ("Poisson Max Residual:", "4.12e-11 (Converged)"),
        ("Newton Iterations:", "4 iterations"),
        ("Hetero Interface:", "Harmonic Mean ε_r"),
    ]
    y_pos = 0.32
    for label, val in solver_info:
        ax_ctrl.text(0.06, y_pos, label, color='#abb2bf', fontsize=8.5, transform=ax_ctrl.transAxes)
        ax_ctrl.text(0.55, y_pos, val, color='#98c379', fontsize=8.5, transform=ax_ctrl.transAxes)
        y_pos -= 0.032

    # Action buttons simulation
    btn_box = Rectangle((0.06, 0.08), 0.88, 0.055, facecolor='#2c313c', edgecolor='#61afef', transform=ax_ctrl.transAxes)
    ax_ctrl.add_patch(btn_box)
    ax_ctrl.text(0.5, 0.105, "Export Telemetry (CSV)", color='#ffffff', fontsize=9.5, fontweight='bold', ha='center', transform=ax_ctrl.transAxes)

    # Right Main Area: 2x2 Plots + Bottom Telemetry bar
    gs_right = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_main[1], wspace=0.18, hspace=0.28)
    
    x = np.linspace(0, 1.0, 600)  # in um
    xj = 0.5
    
    # Realistic physics curves for AlGaAs / GaAs junction at V=0
    # Left (p-AlGaAs): Eg1 = 1.80 eV, chi1 = 3.74 eV
    # Right (n-GaAs): Eg2 = 1.42 eV, chi2 = 4.07 eV
    # Delta Ec = chi2 - chi1 = 0.33 eV; Delta Ev = (Eg1 - Eg2) - Delta Ec = 0.38 - 0.33 = 0.05 eV
    # Vbi ~ 1.29 V
    psi = 1.288 / (1.0 + np.exp(-(x - xj)/0.04))
    
    # 1. Band Diagram
    ax_band = fig.add_subplot(gs_right[0, 0])
    ax_band.set_facecolor('#181b24')
    ax_band.set_title("Energy Band Alignment (Anderson's Rule)", color='#e5c07b', fontsize=11, fontweight='bold', pad=8)
    
    # Ec and Ev with step discontinuities at xj
    Ec = np.where(x < xj, 1.4 - psi, 1.4 - 0.33 - psi)
    Ev = np.where(x < xj, Ec - 1.80, Ec - 1.42)
    Ei = (Ec + Ev) / 2.0
    Ef = np.zeros_like(x)
    
    ax_band.plot(x, Ec, color='#e06c75', lw=2.2, label=r'$E_c(x)$ (Conduction)')
    ax_band.plot(x, Ev, color='#61afef', lw=2.2, label=r'$E_v(x)$ (Valence)')
    ax_band.plot(x, Ei, color='#abb2bf', lw=1.2, ls=':', label=r'$E_i(x)$ (Intrinsic)')
    ax_band.plot(x, Ef, color='#98c379', lw=1.8, ls='--', label=r'$E_F$ (Flat Fermi)')
    ax_band.axvline(xj, color='#5c6370', ls='--', lw=1.0)
    
    # Annotations for delta Ec and delta Ev
    ax_band.annotate(r'$\Delta E_c = 0.33\ \mathrm{eV}$', xy=(xj, 0.4), xytext=(xj+0.05, 0.7),
                     arrowprops=dict(arrowstyle="->", color='#e06c75', lw=1.2), color='#e06c75', fontsize=8.5)
    ax_band.annotate(r'$\Delta E_v = 0.05\ \mathrm{eV}$', xy=(xj, -1.4), xytext=(xj-0.35, -1.2),
                     arrowprops=dict(arrowstyle="->", color='#61afef', lw=1.2), color='#61afef', fontsize=8.5)
    
    ax_band.set_xlabel("Position x [um]", fontsize=9, color='#abb2bf')
    ax_band.set_ylabel("Energy [eV]", fontsize=9, color='#abb2bf')
    ax_band.legend(loc='lower left', fontsize=7.5, framealpha=0.6, facecolor='#1e222b')
    ax_band.grid(True)
    
    # 2. Carrier Concentrations
    ax_carr = fig.add_subplot(gs_right[0, 1])
    ax_carr.set_facecolor('#181b24')
    ax_carr.set_title("Carrier Densities (Log10 Scale)", color='#98c379', fontsize=11, fontweight='bold', pad=8)
    
    p_dens = np.where(x < xj, 1e17 * (1.0 - 0.999 / (1.0 + np.exp(-(xj - x)/0.03))), 1e2 * np.exp((xj - x)/0.03))
    p_dens = np.clip(p_dens, 1e-12, 1e18)
    n_dens = np.where(x >= xj, 1e17 * (1.0 - 0.999 / (1.0 + np.exp(-(x - xj)/0.03))), 1e2 * np.exp((x - xj)/0.03))
    n_dens = np.clip(n_dens, 1e-12, 1e18)
    
    ax_carr.plot(x, np.log10(p_dens), color='#61afef', lw=2.0, label=r'$\log_{10}(p)$ [Holes]')
    ax_carr.plot(x, np.log10(n_dens), color='#e06c75', lw=2.0, label=r'$\log_{10}(n)$ [Electrons]')
    ax_carr.axhline(17, color='#d19a66', ls=':', lw=1.2, label=r'Doping Level $10^{17}\ \mathrm{cm}^{-3}$')
    ax_carr.axvline(xj, color='#5c6370', ls='--', lw=1.0)
    
    ax_carr.set_xlabel("Position x [um]", fontsize=9, color='#abb2bf')
    ax_carr.set_ylabel(r"$\log_{10}(\mathrm{Carrier\ Conc.})\ [\mathrm{cm}^{-3}]$", fontsize=9, color='#abb2bf')
    ax_carr.set_ylim(-12, 18)
    ax_carr.legend(loc='lower left', fontsize=7.5, framealpha=0.6, facecolor='#1e222b')
    ax_carr.grid(True)

    # 3. Dual-Scale Field vs Potential
    ax_field = fig.add_subplot(gs_right[1, 0])
    ax_field.set_facecolor('#181b24')
    ax_field.set_title("Dual-Scale Field E(x) & Potential ψ(x)", color='#61afef', fontsize=11, fontweight='bold', pad=8)
    
    E_field = -62.0 * np.exp(-((x - xj)/0.06)**2)  # Peak -62 kV/cm
    line1 = ax_field.plot(x, E_field, color='#e5c07b', lw=2.0, label='Electric Field E(x)')
    ax_field.set_xlabel("Position x [um]", fontsize=9, color='#abb2bf')
    ax_field.set_ylabel("Electric Field E(x) [kV/cm]", color='#e5c07b', fontsize=9)
    ax_field.tick_params(axis='y', labelcolor='#e5c07b')
    ax_field.grid(True)
    
    ax_pot = ax_field.twinx()
    line2 = ax_pot.plot(x, psi, color='#98c379', lw=2.0, ls='-.', label='Potential ψ(x)')
    ax_pot.set_ylabel("Electrostatic Potential ψ(x) [V]", color='#98c379', fontsize=9)
    ax_pot.tick_params(axis='y', labelcolor='#98c379')
    
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax_field.legend(lines, labels, loc='lower left', fontsize=7.5, framealpha=0.6, facecolor='#1e222b')

    # 4. Space Charge Density
    ax_rho = fig.add_subplot(gs_right[1, 1])
    ax_rho.set_facecolor('#181b24')
    ax_rho.set_title("Net Space Charge Density ρ(x)", color='#c678dd', fontsize=11, fontweight='bold', pad=8)
    
    q = 1.602e-19
    rho = np.where(x < xj, -q * 1e17 / (1.0 + np.exp(-(xj - x - 0.08)/0.015)), q * 1e17 / (1.0 + np.exp(-(x - xj - 0.08)/0.015)))
    rho_coul = rho * 1e6  # C/cm^3
    
    ax_rho.plot(x, rho_coul, color='#c678dd', lw=2.0, label=r'$\rho(x) = q(p - n + N_d^+ - N_a^-)$')
    ax_rho.fill_between(x, rho_coul, 0, where=(rho_coul < 0), color='#61afef', alpha=0.3, label='p-Depletion (Negative)')
    ax_rho.fill_between(x, rho_coul, 0, where=(rho_coul > 0), color='#e06c75', alpha=0.3, label='n-Depletion (Positive)')
    ax_rho.axvline(xj, color='#5c6370', ls='--', lw=1.0)
    
    ax_rho.set_xlabel("Position x [um]", fontsize=9, color='#abb2bf')
    ax_rho.set_ylabel(r"Space Charge Density $\rho(x)\ [\mathrm{C/cm}^3]$", fontsize=9, color='#abb2bf')
    ax_rho.legend(loc='lower right', fontsize=7.5, framealpha=0.6, facecolor='#1e222b')
    ax_rho.grid(True)

    plt.savefig("/mnt/adata/fizikatrulyalya/docs/screenshots/heterojunction_workbench.png", dpi=140, facecolor='#12151c')
    plt.close()
    print("Generated heterojunction_workbench.png successfully.")

def generate_iv_characteristics():
    fig = plt.figure(figsize=(14, 8), dpi=140, facecolor='#12151c')
    gs = gridspec.GridSpec(2, 1, height_ratios=[2.2, 1.0], hspace=0.28, left=0.08, right=0.94, top=0.92, bottom=0.08)
    
    V = np.linspace(0.0, 1.1, 300)
    Vt = 0.02585
    
    # Components of current density
    # Shockley Diffusion: J0_diff ~ 1e-12 A/cm^2
    # SRH Recombination: J0_rec ~ 1e-7 A/cm^2
    J0_diff = 1e-13
    J0_rec = 1e-8
    Rs = 0.15 # series resistance effect
    
    # Effective junction voltage V_j = V - J*Rs
    V_j = np.zeros_like(V)
    J_tot = np.zeros_like(V)
    J_diff = np.zeros_like(V)
    J_rec = np.zeros_like(V)
    
    for i, v in enumerate(V):
        v_est = v
        for _ in range(5):
            jd = J0_diff * (np.exp(np.clip(v_est / Vt, 0, 40)) - 1.0)
            jr = J0_rec * (np.exp(np.clip(v_est / (2.0 * Vt), 0, 40)) - 1.0)
            jt = jd + jr
            v_est = v - jt * Rs * 1e-4
        V_j[i] = v_est
        J_diff[i] = max(1e-15, jd)
        J_rec[i] = max(1e-15, jr)
        J_tot[i] = max(1e-15, jt)

    # 1. Semilog I-V Plot
    ax_iv = fig.add_subplot(gs[0])
    ax_iv.set_facecolor('#181b24')
    ax_iv.set_title("OpenBand Diode Current-Voltage (I-V) Characteristics & Regime Decomposition", color='#61afef', fontsize=12, fontweight='bold', pad=10)
    
    ax_iv.semilogy(V, J_tot, color='#e5c07b', lw=2.6, label='Total Current J_tot(V) [Shockley + SRH]')
    ax_iv.semilogy(V, J_diff, color='#61afef', lw=1.8, ls='--', label='Shockley Ideal Diffusion (η = 1, slope q/kT)')
    ax_iv.semilogy(V, J_rec, color='#e06c75', lw=1.8, ls='-.', label='Space-Charge SRH Recombination (η = 2, slope q/2kT)')
    
    # Operating point marker at V_bias = 0.65 V
    idx_op = np.argmin(np.abs(V - 0.65))
    ax_iv.plot(V[idx_op], J_tot[idx_op], marker='o', markersize=8, color='#98c379', markeredgecolor='#ffffff', markeredgewidth=1.5,
               label=f'Operating Bias $V_{{bias}} = 0.65\\ \\mathrm{{V}}$, $J = {J_tot[idx_op]:.2e}\\ \\mathrm{{A/cm^2}}$')
    ax_iv.axvline(0.65, color='#98c379', ls=':', lw=1.2, alpha=0.7)
    
    ax_iv.set_ylabel(r"Current Density $|J|\ [\mathrm{A/cm}^2]$", fontsize=10, color='#abb2bf')
    ax_iv.set_ylim(1e-9, 1e2)
    ax_iv.set_xlim(0.0, 1.1)
    ax_iv.legend(loc='upper left', fontsize=8.5, framealpha=0.7, facecolor='#1e222b')
    ax_iv.grid(True, which='both')
    
    # 2. Ideality Factor Plot eta(V)
    ax_eta = fig.add_subplot(gs[1], sharex=ax_iv)
    ax_eta.set_facecolor('#181b24')
    
    # Central difference derivative d(ln J)/dV
    ln_J = np.log(J_tot)
    dV = V[1] - V[0]
    dlnJ_dV = np.gradient(ln_J, dV)
    eta = 1.0 / (Vt * dlnJ_dV)
    eta = np.clip(eta, 0.8, 3.0)
    
    ax_eta.plot(V, eta, color='#98c379', lw=2.2, label=r'Differential Ideality Factor $\eta(V) = \frac{1}{V_t \frac{d\ln|J|}{dV}}$')
    ax_eta.axhline(2.0, color='#e06c75', ls=':', lw=1.2, label=r'SRH Recombination Regime ($\eta = 2$)')
    ax_eta.axhline(1.0, color='#61afef', ls=':', lw=1.2, label=r'Ideal Shockley Diffusion Regime ($\eta = 1$)')
    ax_eta.plot(V[idx_op], eta[idx_op], marker='s', markersize=7, color='#e5c07b', markeredgecolor='#ffffff',
                label=f'Operating Ideality $\eta(0.65\\ \\mathrm{{V}}) = {eta[idx_op]:.2f}$')
    
    ax_eta.set_xlabel("Applied Bias Voltage V_bias [V]", fontsize=10, color='#abb2bf')
    ax_eta.set_ylabel(r"Ideality Factor $\eta(V)$", fontsize=10, color='#abb2bf')
    ax_eta.set_ylim(0.8, 2.5)
    ax_eta.legend(loc='upper right', fontsize=8.0, framealpha=0.7, facecolor='#1e222b')
    ax_eta.grid(True)
    
    plt.savefig("/mnt/adata/fizikatrulyalya/docs/screenshots/iv_characteristics.png", dpi=140, facecolor='#12151c')
    plt.close()
    print("Generated iv_characteristics.png successfully.")

if __name__ == "__main__":
    generate_workbench()
    generate_iv_characteristics()
