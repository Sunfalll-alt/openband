//! Analytical vs Numerical metrics extraction for semiconductor p-n junctions and heterojunctions.
//!
//! Provides Depletion Approximation metrics (Vbi, Wdep, Emax, Cj, band offsets)
//! and compares them with exact numerical Finite Difference Poisson-Boltzmann solutions.

use serde::{Deserialize, Serialize};

use super::materials::{EPS_0, MaterialProperties, Q_ELEM, thermal_voltage};

/// Comparison telemetry between analytical depletion model and numerical simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeviceMetrics {
    // --- Analytical (Depletion Approximation) ---
    /// Analytical built-in potential V_bi in Volts.
    pub v_bi_analytical: f64,
    /// Analytical total depletion width W_dep in micrometers (um).
    pub w_dep_analytical_um: f64,
    /// Analytical p-side depletion width x_p in micrometers (um).
    pub xp_analytical_um: f64,
    /// Analytical n-side depletion width x_n in micrometers (um).
    pub xn_analytical_um: f64,
    /// Analytical maximum electric field magnitude at the metallurgical junction in kV/cm.
    pub e_max_analytical_kv_cm: f64,
    /// Analytical junction capacitance per unit area C_j = eps / W_dep in nF/cm^2.
    pub c_j_analytical_nf_cm2: f64,

    // --- Numerical (Finite-Difference Poisson-Boltzmann) ---
    /// Numerical built-in potential in Volts: (psi(L) - psi(0)) + V_bias.
    pub v_bi_numerical: f64,
    /// Numerical potential drop across device contacts psi(L) - psi(0) in Volts.
    pub delta_psi_numerical: f64,
    /// Numerical space-charge depletion width in micrometers (um) extracted from E-field.
    pub w_dep_numerical_um: f64,
    /// Numerical maximum electric field magnitude in kV/cm.
    pub e_max_numerical_kv_cm: f64,
    /// Numerical junction capacitance per unit area in nF/cm^2.
    pub c_j_numerical_nf_cm2: f64,

    // --- Heterojunction Band Offsets ---
    /// Conduction band offset Delta Ec = chi_n - chi_p in eV (at junction).
    pub delta_ec_ev: f64,
    /// Valence band offset Delta Ev = Delta Ec + Eg_n - Eg_p in eV (at junction).
    pub delta_ev_ev: f64,
}

impl DeviceMetrics {
    /// Computes analytical metrics using the abrupt junction depletion approximation.
    pub fn compute_analytical(
        props_p: &MaterialProperties,
        props_n: &MaterialProperties,
        temp_k: f64,
        na: f64,
        nd: f64,
        v_bias: f64,
    ) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
        let vt = thermal_voltage(temp_k);
        let ni_p = props_p.ni(temp_k);
        let ni_n = props_n.ni(temp_k);
        let eps_p = EPS_0 * props_p.eps_r;
        let eps_n = EPS_0 * props_n.eps_r;
        let eps_eff = 2.0 * eps_p * eps_n / (eps_p + eps_n);

        let eg_p = props_p.bandgap(temp_k);
        let eg_n = props_n.bandgap(temp_k);
        let nc_p = props_p.nc(temp_k);
        let nv_p = props_p.nv(temp_k);
        let nc_n = props_n.nc(temp_k);
        let nv_n = props_n.nv(temp_k);

        // Band offsets: Ec(x_j-) - Ec(x_j+) = chi_n - chi_p
        let delta_ec = props_n.electron_affinity - props_p.electron_affinity;
        let delta_ev = delta_ec + eg_n - eg_p;

        // Generalized Anderson heterojunction built-in potential:
        // V_bi = (chi_p - chi_n)/q + (Eg_p - Eg_n)/(2q) + Vt * ln(Na * Nd / (ni_p * ni_n)) + (Vt/2) * ln((Nv_n * Nc_p) / (Nv_p * Nc_n))
        let delta_chi = props_p.electron_affinity - props_n.electron_affinity;
        let delta_eg_half = 0.5 * (eg_p - eg_n);
        let doping_term = vt * ((na * nd) / (ni_p * ni_n)).ln();
        let dos_term = 0.5 * vt * ((nv_n * nc_p) / (nv_p * nc_n)).ln();
        let v_bi = (delta_chi + delta_eg_half + doping_term + dos_term).max(0.0);

        // Effective barrier under bias
        let v_barrier = (v_bi - v_bias).max(0.001);

        // Heterojunction depletion width:
        // W_dep = sqrt( (2 * eps_p * eps_n * (Na + Nd)^2 * V_barrier) / (q * Na * Nd * (eps_p * Na + eps_n * Nd)) )
        // Or in effective epsilon:
        let w_dep_cm = ((2.0 * eps_eff * v_barrier / Q_ELEM) * (1.0 / na + 1.0 / nd)).sqrt();
        let w_dep_um = w_dep_cm * 1e4;

        let xp_um = w_dep_um * (nd / (na + nd));
        let xn_um = w_dep_um * (na / (na + nd));

        // Peak field |E_max| = 2 * V_barrier / W_dep in V/cm -> convert to kV/cm
        let e_max_v_cm = 2.0 * v_barrier / w_dep_cm.max(1e-12);
        let e_max_kv_cm = e_max_v_cm * 1e-3;

        // C_j = eps_eff / W_dep in F/cm^2 -> convert to nF/cm^2
        let c_j_f_cm2 = eps_eff / w_dep_cm.max(1e-12);
        let c_j_nf_cm2 = c_j_f_cm2 * 1e9;

        (
            v_bi,
            w_dep_um,
            xp_um,
            xn_um,
            e_max_kv_cm,
            c_j_nf_cm2,
            delta_ec,
            delta_ev,
        )
    }

    /// Evaluates both analytical and numerical metrics from simulation solution arrays.
    pub fn extract(
        props_p: &MaterialProperties,
        props_n: &MaterialProperties,
        temp_k: f64,
        na: f64,
        nd: f64,
        v_bias: f64,
        x_cm: &[f64],
        psi_v: &[f64],
        e_field_kv_cm: &[f64],
    ) -> Self {
        let (
            v_bi_ana,
            w_dep_ana_um,
            xp_ana_um,
            xn_ana_um,
            e_max_ana_kv_cm,
            c_j_ana_nf_cm2,
            delta_ec,
            delta_ev,
        ) = Self::compute_analytical(props_p, props_n, temp_k, na, nd, v_bias);

        let delta_psi = if !psi_v.is_empty() {
            psi_v[psi_v.len() - 1] - psi_v[0]
        } else {
            0.0
        };
        let v_bi_num = delta_psi + v_bias;

        // Peak electric field magnitude
        let mut e_max_num_kv_cm = 0.0;
        for &e in e_field_kv_cm {
            let e_abs = e.abs();
            if e_abs > e_max_num_kv_cm {
                e_max_num_kv_cm = e_abs;
            }
        }

        // Numerical depletion width extraction:
        // Find interval [x_start, x_end] where |E(x)| >= 2% of |E_max|
        let threshold = 0.02 * e_max_num_kv_cm;
        let mut first_idx = None;
        let mut last_idx = None;

        for (i, &e) in e_field_kv_cm.iter().enumerate() {
            if e.abs() >= threshold {
                if first_idx.is_none() {
                    first_idx = Some(i);
                }
                last_idx = Some(i);
            }
        }

        let w_dep_num_um = match (first_idx, last_idx) {
            (Some(i1), Some(i2)) if i2 >= i1 => {
                let dx_cm = x_cm[i2] - x_cm[i1];
                dx_cm * 1e4
            }
            _ => 0.0,
        };

        let eps_eff = 2.0 * EPS_0 * props_p.eps_r * props_n.eps_r / (props_p.eps_r + props_n.eps_r);
        let w_dep_num_cm = w_dep_num_um * 1e-4;
        let c_j_num_nf_cm2 = if w_dep_num_cm > 1e-12 {
            (eps_eff / w_dep_num_cm) * 1e9
        } else {
            0.0
        };

        Self {
            v_bi_analytical: v_bi_ana,
            w_dep_analytical_um: w_dep_ana_um,
            xp_analytical_um: xp_ana_um,
            xn_analytical_um: xn_ana_um,
            e_max_analytical_kv_cm: e_max_ana_kv_cm,
            c_j_analytical_nf_cm2: c_j_ana_nf_cm2,
            v_bi_numerical: v_bi_num,
            delta_psi_numerical: delta_psi,
            w_dep_numerical_um: w_dep_num_um,
            e_max_numerical_kv_cm: e_max_num_kv_cm,
            c_j_numerical_nf_cm2: c_j_num_nf_cm2,
            delta_ec_ev: delta_ec,
            delta_ev_ev: delta_ev,
        }
    }
}
