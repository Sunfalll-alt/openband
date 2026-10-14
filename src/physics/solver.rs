//! 1D Non-Linear Poisson-Boltzmann Solver using Finite-Difference Method (FDM),
//! Damped Newton-Raphson, and Tridiagonal Matrix Algorithm (TDMA / Thomas algorithm).
//!
//! Solves:
//!     d/dx [ eps_r(x) * dpsi/dx ] = - q / eps_0 * [ p(x) - n(x) + Nd+(x) - Na-(x) ]
//! Supporting both homojunctions and heterojunctions with Anderson's electron affinity rule.

use serde::{Deserialize, Serialize};

use super::iv_curve::{IvCurveResult, compute_iv_curve};
use super::materials::{MaterialProperties, MaterialType, Q_ELEM, thermal_voltage};
use super::metrics::DeviceMetrics;

/// Standard spatial grid node count.
pub const NUM_NODES: usize = 600;

/// Device configuration parameters for 1D simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeviceParams {
    /// In homojunction mode: material for entire device.
    /// In heterojunction mode: p-side material (x < xj).
    pub material: MaterialType,
    /// n-side material (x >= xj) when is_heterojunction is true.
    pub material_n: MaterialType,
    /// Whether heterojunction mode (different p and n materials) is enabled.
    pub is_heterojunction: bool,
    /// Absolute temperature in Kelvin (100 K to 600 K).
    pub temp_k: f64,
    /// Acceptor doping concentration Na in cm^-3 (1e14 to 1e19).
    pub na: f64,
    /// Donor doping concentration Nd in cm^-3 (1e14 to 1e19).
    pub nd: f64,
    /// Total device length in micrometers (0.2 um to 5.0 um).
    pub length_um: f64,
    /// Metallurgical junction position as fraction of length (0.1 to 0.9).
    pub junction_ratio: f64,
    /// Applied external bias in Volts (-5.0 V to +0.8 V).
    pub v_bias: f64,
    /// Carrier lifetime tau_0 in seconds (for SRH and diffusion, 1e-9 to 1e-5 s).
    pub tau_0: f64,
}

impl Default for DeviceParams {
    fn default() -> Self {
        Self {
            material: MaterialType::Silicon,
            material_n: MaterialType::Silicon,
            is_heterojunction: false,
            temp_k: 300.0,
            na: 1.0e16,
            nd: 1.0e16,
            length_um: 1.0,
            junction_ratio: 0.5,
            v_bias: 0.0,
            tau_0: 1.0e-7, // 100 ns
        }
    }
}

/// Simulation solution containing spatial distribution profiles and comparison metrics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationResult {
    /// Spatial grid coordinates in micrometers (um).
    pub x_um: Vec<f64>,
    /// Electrostatic potential psi(x) in Volts.
    pub psi_v: Vec<f64>,
    /// Electric field E(x) = -dpsi/dx in kV/cm.
    pub e_field_kv_cm: Vec<f64>,
    /// Net space charge density rho(x) in C/cm^3.
    pub rho_c_cm3: Vec<f64>,
    /// Free electron concentration n(x) in cm^-3.
    pub n_cm3: Vec<f64>,
    /// Free hole concentration p(x) in cm^-3.
    pub p_cm3: Vec<f64>,
    /// Conduction band edge Ec(x) in eV.
    pub e_c_ev: Vec<f64>,
    /// Valence band edge Ev(x) in eV.
    pub e_v_ev: Vec<f64>,
    /// Intrinsic Fermi level Ei(x) in eV.
    pub e_i_ev: Vec<f64>,
    /// Electron quasi-Fermi level E_Fn(x) in eV.
    pub e_fn_ev: Vec<f64>,
    /// Hole quasi-Fermi level E_Fp(x) in eV.
    pub e_fp_ev: Vec<f64>,
    /// Net doping profile (Nd - Na) in cm^-3.
    pub net_doping_cm3: Vec<f64>,
    /// Iterations taken by Newton-Raphson.
    pub iterations: usize,
    /// Convergence flag.
    pub converged: bool,
    /// Final maximum potential step update max |delta_u|.
    pub max_delta_u: f64,
    /// Extracted device metrics (analytical vs numerical).
    pub metrics: DeviceMetrics,
    /// Computed I-V characteristic curve.
    pub iv_curve: IvCurveResult,
}

/// Zero-allocation scratch workspace for Newton-Raphson solver iterations.
pub struct SolverWorkspace {
    pub x_cm: Vec<f64>,
    pub u: Vec<f64>,
    pub phi_fn: Vec<f64>,
    pub phi_fp: Vec<f64>,
    pub c_doping: Vec<f64>,
    pub kappa: Vec<f64>,
    pub theta_shift: Vec<f64>,
    pub gamma_ni: Vec<f64>,
    pub residual: Vec<f64>,
    pub diag_a: Vec<f64>,
    pub diag_b: Vec<f64>,
    pub diag_c: Vec<f64>,
    pub tdma_c_prime: Vec<f64>,
    pub tdma_d_prime: Vec<f64>,
    pub delta_u: Vec<f64>,
}

impl Default for SolverWorkspace {
    fn default() -> Self {
        Self::new(NUM_NODES)
    }
}

impl SolverWorkspace {
    /// Allocates reusable scratch buffers for N grid points.
    pub fn new(n: usize) -> Self {
        let m = if n >= 2 { n - 2 } else { 0 };
        Self {
            x_cm: vec![0.0; n],
            u: vec![0.0; n],
            phi_fn: vec![0.0; n],
            phi_fp: vec![0.0; n],
            c_doping: vec![0.0; n],
            kappa: vec![1.0; n],
            theta_shift: vec![0.0; n],
            gamma_ni: vec![1.0; n],
            residual: vec![0.0; m],
            diag_a: vec![0.0; m],
            diag_b: vec![0.0; m],
            diag_c: vec![0.0; m],
            tdma_c_prime: vec![0.0; m],
            tdma_d_prime: vec![0.0; m],
            delta_u: vec![0.0; m],
        }
    }

    /// Solves the 1D Poisson-Boltzmann equation for the given device parameters.
    ///
    /// Performs ZERO heap allocations during the Newton-Raphson iterations.
    pub fn solve(&mut self, params: &DeviceParams) -> SimulationResult {
        let n = NUM_NODES;
        let m = n - 2;

        let props_p: MaterialProperties = params.material.properties();
        let props_n: MaterialProperties = if params.is_heterojunction {
            params.material_n.properties()
        } else {
            props_p
        };

        let temp_k = params.temp_k.max(10.0);
        let vt = thermal_voltage(temp_k);
        let ni_p = props_p.ni(temp_k);
        let ni_n = props_n.ni(temp_k);
        let na = params.na.max(1.0);
        let nd = params.nd.max(1.0);
        let length_cm = params.length_um.max(0.01) * 1e-4;
        let dx_cm = length_cm / ((n - 1) as f64);
        let xj_cm = (params.junction_ratio.clamp(0.01, 0.99)) * length_cm;

        // Spatial grid setup
        for i in 0..n {
            self.x_cm[i] = (i as f64) * dx_cm;
        }

        // Reference scaling: p-side material at x=0
        let debye_l_cm = props_p.debye_length(temp_k);
        let h = dx_cm / debye_l_cm;
        let h2 = h * h;
        let inv_h2 = 1.0 / h2;

        // Band alignment shift Delta theta_n across heterojunction:
        // theta_n = [ (chi_n - chi_p) + 0.5*(Eg_n - Eg_p) + 0.5*Vt*ln(Nc_n*Nv_p / (Nv_n*Nc_p)) ] / Vt
        let eg_p = props_p.bandgap(temp_k);
        let eg_n = props_n.bandgap(temp_k);
        let nc_p = props_p.nc(temp_k);
        let nv_p = props_p.nv(temp_k);
        let nc_n = props_n.nc(temp_k);
        let nv_n = props_n.nv(temp_k);

        let delta_chi = props_n.electron_affinity - props_p.electron_affinity;
        let delta_eg = eg_n - eg_p;
        let dos_ratio_term = 0.5 * vt * ((nc_n * nv_p) / (nv_n * nc_p)).ln();
        let delta_theta_n = if params.is_heterojunction {
            (delta_chi + 0.5 * delta_eg + dos_ratio_term) / vt
        } else {
            0.0
        };

        // Material distributions across the grid
        for i in 0..n {
            let x = self.x_cm[i];
            if x < xj_cm {
                self.kappa[i] = 1.0; // eps_r / eps_r_p
                self.theta_shift[i] = 0.0;
                self.gamma_ni[i] = 1.0;
                self.c_doping[i] = -na / ni_p;
            } else {
                self.kappa[i] = props_n.eps_r / props_p.eps_r;
                self.theta_shift[i] = delta_theta_n;
                self.gamma_ni[i] = ni_n / ni_p;
                self.c_doping[i] = nd / ni_p;
            }
        }

        // Boundary conditions at Ohmic contacts
        let u_0 = -((na / (2.0 * ni_p)).asinh());
        let u_end = ((nd / (2.0 * ni_n)).asinh() - delta_theta_n) - (params.v_bias / vt);

        // Smoothly interpolated quasi-Fermi potentials across metallurgical junction
        let v_bias_norm = -params.v_bias / vt;
        for i in 0..n {
            let x = self.x_cm[i];
            if x <= xj_cm {
                let s = (x / xj_cm).clamp(0.0, 1.0);
                let s_step = s * s * (3.0 - 2.0 * s);
                self.phi_fn[i] = v_bias_norm * s_step;
                self.phi_fp[i] = 0.0;
            } else {
                let tau = ((x - xj_cm) / (length_cm - xj_cm)).clamp(0.0, 1.0);
                let tau_step = tau * tau * (3.0 - 2.0 * tau);
                self.phi_fn[i] = v_bias_norm;
                self.phi_fp[i] = v_bias_norm * tau_step;
            }
        }

        // Smooth initial guess connecting boundary conditions u_0 and u_end
        let w_init_cm = (0.05 * length_cm).max(2.0 * debye_l_cm);
        for i in 0..n {
            let x = self.x_cm[i];
            let arg = -((x - xj_cm) / w_init_cm).clamp(-40.0, 40.0);
            let s = 1.0 / (1.0 + arg.exp());
            self.u[i] = u_0 + (u_end - u_0) * s;
        }
        self.u[0] = u_0;
        self.u[n - 1] = u_end;

        // Newton-Raphson Solver Loop
        let max_iters = 100;
        let tolerance = 1e-8;
        let mut converged = false;
        let mut final_max_delta_u = 0.0;
        let mut iterations_done = 0;

        for iter in 0..max_iters {
            iterations_done = iter + 1;

            // 1. Build Residual F_i and Tridiagonal Jacobian J
            for k in 0..m {
                let i = k + 1;
                let u_i = self.u[i];
                let u_prev = self.u[i - 1];
                let u_next = self.u[i + 1];

                let u_eff = u_i + self.theta_shift[i];
                let phi_fn_i = self.phi_fn[i];
                let phi_fp_i = self.phi_fp[i];
                let c_i = self.c_doping[i];
                let gamma_i = self.gamma_ni[i];

                let exp_n_arg = (u_eff - phi_fn_i).clamp(-600.0, 600.0);
                let exp_p_arg = (phi_fp_i - u_eff).clamp(-600.0, 600.0);

                let n_norm = gamma_i * exp_n_arg.exp();
                let p_norm = gamma_i * exp_p_arg.exp();

                // Spatially varying permittivity:
                // kappa_{i-1/2} = 2 * kappa_{i-1} * kappa_i / (kappa_{i-1} + kappa_i)
                // kappa_{i+1/2} = 2 * kappa_i * kappa_{i+1} / (kappa_i + kappa_{i+1})
                let k_prev = self.kappa[i - 1];
                let k_cur = self.kappa[i];
                let k_next = self.kappa[i + 1];

                let k_minus_half = (2.0 * k_prev * k_cur) / (k_prev + k_cur);
                let k_plus_half = (2.0 * k_cur * k_next) / (k_cur + k_next);

                // Lap = [ kappa_{i-1/2} * u_{i-1} - (kappa_{i-1/2} + kappa_{i+1/2}) * u_i + kappa_{i+1/2} * u_{i+1} ] / h^2
                let laplacian = (k_minus_half * u_prev - (k_minus_half + k_plus_half) * u_i
                    + k_plus_half * u_next)
                    * inv_h2;

                let residual_i = laplacian - (n_norm - p_norm - c_i);

                self.residual[k] = residual_i;
                self.diag_a[k] = k_minus_half * inv_h2;
                self.diag_c[k] = k_plus_half * inv_h2;
                self.diag_b[k] = -(k_minus_half + k_plus_half) * inv_h2 - (n_norm + p_norm);
            }

            // 2. TDMA (Thomas Algorithm) to solve J * delta_u = -F
            // Forward sweep
            let b0 = self.diag_b[0];
            self.tdma_c_prime[0] = self.diag_c[0] / b0;
            self.tdma_d_prime[0] = (-self.residual[0]) / b0;

            for k in 1..m {
                let a_k = self.diag_a[k];
                let denom = self.diag_b[k] - a_k * self.tdma_c_prime[k - 1];
                self.tdma_c_prime[k] = self.diag_c[k] / denom;
                self.tdma_d_prime[k] = (-self.residual[k] - a_k * self.tdma_d_prime[k - 1]) / denom;
            }

            // Back substitution
            self.delta_u[m - 1] = self.tdma_d_prime[m - 1];
            for k in (0..m - 1).rev() {
                self.delta_u[k] = self.tdma_d_prime[k] - self.tdma_c_prime[k] * self.delta_u[k + 1];
            }

            // 3. Convergence check & Damped updates (Bank-Rose / step clamping <= 2.0)
            let mut max_step: f64 = 0.0;
            for k in 0..m {
                let step_abs = self.delta_u[k].abs();
                if step_abs > max_step {
                    max_step = step_abs;
                }
            }

            final_max_delta_u = max_step;

            if max_step < tolerance {
                converged = true;
                break;
            }

            for k in 0..m {
                let i = k + 1;
                let step = self.delta_u[k].clamp(-2.0, 2.0);
                self.u[i] += step;
            }
        }

        // Ensure boundary conditions remain strictly enforced
        self.u[0] = u_0;
        self.u[n - 1] = u_end;

        // Post-processing to physical units
        let mut x_um = Vec::with_capacity(n);
        let mut psi_v = Vec::with_capacity(n);
        let mut e_field_kv_cm = Vec::with_capacity(n);
        let mut rho_c_cm3 = Vec::with_capacity(n);
        let mut n_cm3 = Vec::with_capacity(n);
        let mut p_cm3 = Vec::with_capacity(n);
        let mut e_c_ev = Vec::with_capacity(n);
        let mut e_v_ev = Vec::with_capacity(n);
        let mut e_i_ev = Vec::with_capacity(n);
        let mut e_fn_ev = Vec::with_capacity(n);
        let mut e_fp_ev = Vec::with_capacity(n);
        let mut net_doping_cm3 = Vec::with_capacity(n);

        // Anderson's vacuum level reference constant:
        // Set Evac so that at x=0 (p-neutral region): Ec(0) = -psi(0) + 0.5*Eg_p + 0.5*Vt*ln(Nc_p/Nv_p)
        let evac_ref = props_p.electron_affinity + 0.5 * eg_p + 0.5 * vt * (nc_p / nv_p).ln();

        for i in 0..n {
            let x_c = self.x_cm[i];
            x_um.push(x_c * 1e4);

            let psi = self.u[i] * vt;
            psi_v.push(psi);

            let cur_props = if x_c < xj_cm { &props_p } else { &props_n };
            let cur_ni = cur_props.ni(temp_k);
            let cur_eg = cur_props.bandgap(temp_k);
            let cur_nc = cur_props.nc(temp_k);
            let cur_nv = cur_props.nv(temp_k);
            let cur_chi = cur_props.electron_affinity;

            let u_eff = self.u[i] + self.theta_shift[i];
            let exp_n = (u_eff - self.phi_fn[i]).clamp(-600.0, 600.0).exp();
            let exp_p = (self.phi_fp[i] - u_eff).clamp(-600.0, 600.0).exp();

            let n_val = (cur_ni * exp_n).max(1e-20);
            let p_val = (cur_ni * exp_p).max(1e-20);
            n_cm3.push(n_val);
            p_cm3.push(p_val);

            let c_val = self.c_doping[i] * ni_p;
            net_doping_cm3.push(c_val);

            // rho(x) = q * [ p(x) - n(x) + Nd+(x) - Na-(x) ] = q * [ p(x) - n(x) + c_val ]
            let rho = Q_ELEM * (p_val - n_val + c_val);
            rho_c_cm3.push(rho);

            // Anderson's rule for energy bands:
            // Evac(x) = - psi(x) + evac_ref
            // Ec(x) = Evac(x) - chi(x)
            // Ev(x) = Ec(x) - Eg(x)
            let e_vac = -psi + evac_ref;
            let ec = e_vac - cur_chi;
            e_c_ev.push(ec);

            let ev = ec - cur_eg;
            e_v_ev.push(ev);

            let ei = ec - 0.5 * cur_eg - 0.5 * vt * (cur_nc / cur_nv).ln();
            e_i_ev.push(ei);

            // Quasi-Fermi levels in eV
            e_fn_ev.push(-self.phi_fn[i] * vt);
            e_fp_ev.push(-self.phi_fp[i] * vt);
        }

        // Electric field E(x) = - dpsi / dx via finite differences (kV/cm)
        for i in 0..n {
            let dpsi_dx = if i == 0 {
                (psi_v[1] - psi_v[0]) / dx_cm
            } else if i == n - 1 {
                (psi_v[n - 1] - psi_v[n - 2]) / dx_cm
            } else {
                (psi_v[i + 1] - psi_v[i - 1]) / (2.0 * dx_cm)
            };
            let e_field_v_cm = -dpsi_dx;
            e_field_kv_cm.push(e_field_v_cm * 1e-3);
        }

        // Extract comparison telemetry
        let metrics = DeviceMetrics::extract(
            &props_p,
            &props_n,
            temp_k,
            na,
            nd,
            params.v_bias,
            &self.x_cm,
            &psi_v,
            &e_field_kv_cm,
        );

        // Compute I-V curve
        let iv_curve = compute_iv_curve(
            &props_p,
            &props_n,
            temp_k,
            na,
            nd,
            params.v_bias,
            params.tau_0,
            metrics.v_bi_numerical,
        );

        SimulationResult {
            x_um,
            psi_v,
            e_field_kv_cm,
            rho_c_cm3,
            n_cm3,
            p_cm3,
            e_c_ev,
            e_v_ev,
            e_i_ev,
            e_fn_ev,
            e_fp_ev,
            net_doping_cm3,
            iterations: iterations_done,
            converged,
            max_delta_u: final_max_delta_u,
            metrics,
            iv_curve,
        }
    }
}

/// Convenience standalone function to solve Poisson-Boltzmann.
pub fn solve(params: &DeviceParams) -> SimulationResult {
    let mut ws = SolverWorkspace::new(NUM_NODES);
    ws.solve(params)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_silicon_equilibrium_convergence() {
        let params = DeviceParams {
            material: MaterialType::Silicon,
            material_n: MaterialType::Silicon,
            is_heterojunction: false,
            temp_k: 300.0,
            na: 1.0e16,
            nd: 1.0e16,
            length_um: 1.0,
            junction_ratio: 0.5,
            v_bias: 0.0,
            tau_0: 1.0e-7,
        };

        let result = solve(&params);

        assert!(
            result.converged,
            "Solver failed to converge! Final max_delta_u = {}",
            result.max_delta_u
        );
        assert!(
            result.iterations < 50,
            "Iterations ({}) exceeded target limit of 50",
            result.iterations
        );
    }

    #[test]
    fn test_silicon_built_in_potential_analytical_match() {
        let params = DeviceParams {
            material: MaterialType::Silicon,
            material_n: MaterialType::Silicon,
            is_heterojunction: false,
            temp_k: 300.0,
            na: 1.0e16,
            nd: 1.0e16,
            length_um: 1.0,
            junction_ratio: 0.5,
            v_bias: 0.0,
            tau_0: 1.0e-7,
        };

        let result = solve(&params);
        let delta_psi = result.psi_v[NUM_NODES - 1] - result.psi_v[0];
        let v_bi_analytical = result.metrics.v_bi_analytical;

        let relative_error = ((delta_psi - v_bi_analytical).abs()) / v_bi_analytical;
        assert!(
            relative_error < 0.01,
            "Built-in potential delta_psi ({:.4} V) deviated from analytical ({:.4} V) by {:.2}% (> 1%)",
            delta_psi,
            v_bi_analytical,
            relative_error * 100.0
        );
    }

    #[test]
    fn test_charge_neutrality_at_boundaries() {
        let params = DeviceParams {
            material: MaterialType::Silicon,
            material_n: MaterialType::Silicon,
            is_heterojunction: false,
            temp_k: 300.0,
            na: 1.0e16,
            nd: 1.0e16,
            length_um: 1.0,
            junction_ratio: 0.5,
            v_bias: 0.0,
            tau_0: 1.0e-7,
        };

        let result = solve(&params);
        let rho_0 = result.rho_c_cm3[0];
        let rho_end = result.rho_c_cm3[NUM_NODES - 1];

        assert!(
            rho_0.abs() < 1e-6,
            "Charge neutrality violated at p-boundary: rho(0) = {} C/cm^3",
            rho_0
        );
        assert!(
            rho_end.abs() < 1e-6,
            "Charge neutrality violated at n-boundary: rho(L) = {} C/cm^3",
            rho_end
        );
    }

    #[test]
    fn test_heterojunction_algaas_gaas() {
        let params = DeviceParams {
            material: MaterialType::AlGaAs,
            material_n: MaterialType::GalliumArsenide,
            is_heterojunction: true,
            temp_k: 300.0,
            na: 5.0e17,
            nd: 1.0e16,
            length_um: 1.0,
            junction_ratio: 0.5,
            v_bias: 0.0,
            tau_0: 1.0e-7,
        };

        let result = solve(&params);
        assert!(result.converged, "Heterojunction solver failed to converge");
        assert!(
            result.metrics.delta_ec_ev.abs() > 0.01,
            "Delta Ec must be non-zero"
        );
        assert!(
            result.metrics.delta_ev_ev.abs() > 0.01,
            "Delta Ev must be non-zero"
        );

        // Chi(GaAs) = 4.07, Chi(AlGaAs) = 3.74 => Delta Ec = 0.33 eV
        assert!(
            (result.metrics.delta_ec_ev - 0.33).abs() < 0.02,
            "Delta Ec expected ~0.33 eV, got {}",
            result.metrics.delta_ec_ev
        );

        // Analytical vs. Numerical V_bi must match within 1-2%
        let v_bi_ana = result.metrics.v_bi_analytical;
        let v_bi_num = result.metrics.v_bi_numerical;
        let diff_pct = ((v_bi_num - v_bi_ana).abs() / v_bi_ana) * 100.0;
        assert!(
            diff_pct < 2.0,
            "Heterojunction V_bi mismatch: Ana = {:.4} V, Num = {:.4} V (diff = {:.2}%)",
            v_bi_ana,
            v_bi_num,
            diff_pct
        );
    }

    #[test]
    fn test_asymmetric_junction_and_bias() {
        let mut ws = SolverWorkspace::new(NUM_NODES);

        // Reverse bias
        let params_rev = DeviceParams {
            material: MaterialType::Silicon,
            material_n: MaterialType::Silicon,
            is_heterojunction: false,
            temp_k: 300.0,
            na: 1.0e17,
            nd: 1.0e15,
            length_um: 2.0,
            junction_ratio: 0.4,
            v_bias: -2.0,
            tau_0: 1.0e-7,
        };
        let res_rev = ws.solve(&params_rev);
        assert!(
            res_rev.converged,
            "Asymmetric reverse bias failed to converge"
        );

        // Forward bias
        let params_fwd = DeviceParams {
            material: MaterialType::Silicon,
            material_n: MaterialType::Silicon,
            is_heterojunction: false,
            temp_k: 300.0,
            na: 1.0e17,
            nd: 1.0e15,
            length_um: 2.0,
            junction_ratio: 0.4,
            v_bias: 0.5,
            tau_0: 1.0e-7,
        };
        let res_fwd = ws.solve(&params_fwd);
        assert!(
            res_fwd.converged,
            "Asymmetric forward bias failed to converge"
        );
    }

    #[test]
    fn test_all_materials_solve() {
        let mut ws = SolverWorkspace::new(NUM_NODES);
        for mat in MaterialType::ALL {
            let params = DeviceParams {
                material: mat,
                material_n: mat,
                is_heterojunction: false,
                temp_k: 300.0,
                na: 1.0e16,
                nd: 1.0e16,
                length_um: 1.0,
                junction_ratio: 0.5,
                v_bias: 0.0,
                tau_0: 1.0e-7,
            };
            let res = ws.solve(&params);
            assert!(
                res.converged,
                "Material {:?} failed to converge at 300K, iters={}",
                mat, res.iterations
            );
        }
    }
}
