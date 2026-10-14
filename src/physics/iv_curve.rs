//! Diode Current-Voltage (I-V) Characteristic Curve Analysis.
//!
//! Models:
//! 1. Ideal Minority Carrier Diffusion Current (Shockley diode equation):
//!    J_diff(V) = J_0,diff * (exp(q*V / (k_B*T)) - 1)
//! 2. Depletion Region Shockley-Read-Hall (SRH) Recombination / Generation Current:
//!    J_SRH(V) = [q * n_i,eff * W_dep(V) / (2 * tau_0)] * (exp(q*V / (2 * k_B*T)) - 1)
//! 3. High-level injection and series resistance (R_s) roll-off.
//! 4. Extraction of saturation current J0, operating point J(V_bias), and ideality factor eta.

use serde::{Deserialize, Serialize};

use super::materials::{EPS_0, MaterialProperties, Q_ELEM, thermal_voltage};

/// Number of voltage evaluation points in the I-V sweep.
pub const NUM_IV_POINTS: usize = 120;

/// Output dataset from the I-V characteristic sweep.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IvCurveResult {
    /// Applied voltage sweep array V in Volts [-2.0 V, +1.0 V].
    pub v_v: Vec<f64>,
    /// Total current density magnitude |J_total| in A/cm^2.
    pub j_total_a_cm2: Vec<f64>,
    /// Signed total current density in A/cm^2.
    pub j_total_signed_a_cm2: Vec<f64>,
    /// Minority carrier diffusion current component |J_diff| in A/cm^2.
    pub j_diff_a_cm2: Vec<f64>,
    /// SRH recombination / generation component |J_SRH| in A/cm^2.
    pub j_srh_a_cm2: Vec<f64>,
    /// Local ideality factor eta(V).
    pub ideality: Vec<f64>,
    /// Extracted reverse saturation current density J0 in A/cm^2.
    pub j0_saturation: f64,
    /// Operating point current density J(V_bias) in A/cm^2.
    pub j_operating: f64,
    /// Operating point local ideality factor eta(V_bias).
    pub ideality_operating: f64,
}

/// Computes the complete I-V characteristic curve across V in [-2.0 V, +1.0 V].
pub fn compute_iv_curve(
    props_p: &MaterialProperties,
    props_n: &MaterialProperties,
    temp_k: f64,
    na: f64,
    nd: f64,
    v_bias: f64,
    tau_0: f64,
    v_bi: f64,
) -> IvCurveResult {
    let vt = thermal_voltage(temp_k);
    let ni_p = props_p.ni(temp_k);
    let ni_n = props_n.ni(temp_k);
    let ni_eff = (ni_p * ni_n).sqrt();

    let eps_p = EPS_0 * props_p.eps_r;
    let eps_n = EPS_0 * props_n.eps_r;
    let eps_eff = 2.0 * eps_p * eps_n / (eps_p + eps_n);

    // Einstein relation: D = mu * Vt
    let dn = props_p.mu_n * vt; // Minority electron diffusion on p-side
    let dp = props_n.mu_p * vt; // Minority hole diffusion on n-side

    let ln = (dn * tau_0.max(1e-12)).sqrt();
    let lp = (dp * tau_0.max(1e-12)).sqrt();

    // Shockley diffusion saturation current:
    // J_0,diff = q * (Dn/Ln * ni_p^2 / Na + Dp/Lp * ni_n^2 / Nd)
    let j0_diff = Q_ELEM * ((dn / ln) * (ni_p * ni_p / na) + (dp / lp) * (ni_n * ni_n / nd));

    // Lumped series resistance in Ohm * cm^2
    let r_s = 0.02;

    let v_min = -2.0;
    let v_max = 1.0;
    let n_pts = NUM_IV_POINTS;
    let dv = (v_max - v_min) / ((n_pts - 1) as f64);

    let mut v_v = Vec::with_capacity(n_pts);
    let mut j_total_a_cm2 = Vec::with_capacity(n_pts);
    let mut j_total_signed_a_cm2 = Vec::with_capacity(n_pts);
    let mut j_diff_a_cm2 = Vec::with_capacity(n_pts);
    let mut j_srh_a_cm2 = Vec::with_capacity(n_pts);

    // Helper closure to evaluate ideal currents at internal junction voltage Vj
    let eval_currents = |vj: f64| -> (f64, f64) {
        let v_bar = (v_bi - vj).max(0.01);
        let w_dep_cm = ((2.0 * eps_eff * v_bar / Q_ELEM) * (1.0 / na + 1.0 / nd)).sqrt();

        // J_diff
        let exp_diff_arg = (vj / vt).clamp(-40.0, 40.0);
        let j_diff = j0_diff * (exp_diff_arg.exp() - 1.0);

        // J_SRH
        let exp_srh_arg = (vj / (2.0 * vt)).clamp(-40.0, 40.0);
        let j0_srh = (Q_ELEM * ni_eff * w_dep_cm) / (2.0 * tau_0.max(1e-12));
        let j_srh = j0_srh * (exp_srh_arg.exp() - 1.0);

        (j_diff, j_srh)
    };

    // Helper closure to solve for internal junction voltage Vj taking Rs into account:
    // Vj + (J_diff(Vj) + J_SRH(Vj)) * Rs = V_ext
    let solve_vj = |v_ext: f64| -> f64 {
        if v_ext <= 0.0 {
            return v_ext;
        }
        let mut vj = v_ext.min(v_bi + 0.1);
        for _ in 0..6 {
            let (jd, js) = eval_currents(vj);
            let j_tot = jd + js;
            let f = vj + j_tot * r_s - v_ext;
            // Derivative df/dvj:
            let djd_dvj = (j0_diff / vt) * (vj / vt).clamp(-40.0, 40.0).exp();
            let v_bar = (v_bi - vj).max(0.01);
            let w_dep_cm = ((2.0 * eps_eff * v_bar / Q_ELEM) * (1.0 / na + 1.0 / nd)).sqrt();
            let j0_srh = (Q_ELEM * ni_eff * w_dep_cm) / (2.0 * tau_0.max(1e-12));
            let djs_dvj = (j0_srh / (2.0 * vt)) * (vj / (2.0 * vt)).clamp(-40.0, 40.0).exp();
            let df = 1.0 + (djd_dvj + djs_dvj) * r_s;
            let step = (f / df).clamp(-0.2, 0.2);
            vj -= step;
            if step.abs() < 1e-6 {
                break;
            }
        }
        vj
    };

    for i in 0..n_pts {
        let v = v_min + (i as f64) * dv;
        let vj = solve_vj(v);
        let (jd, js) = eval_currents(vj);
        let j_tot = jd + js;

        v_v.push(v);
        j_total_signed_a_cm2.push(j_tot);
        j_total_a_cm2.push(j_tot.abs().max(1e-25));
        j_diff_a_cm2.push(jd.abs().max(1e-25));
        j_srh_a_cm2.push(js.abs().max(1e-25));
    }

    // High-accuracy smoothed ideality factor eta = 1 / (Vt * d(ln J)/dV)
    // Uses central differences for interior points, second-order one-sided differences at boundaries
    let mut ideality = vec![2.0; n_pts];
    let mut ln_j = vec![0.0; n_pts];
    for i in 0..n_pts {
        ln_j[i] = j_total_signed_a_cm2[i].max(1e-25).ln();
    }

    // First index where forward bias is established (V >= 0.05 V)
    let start_idx = v_v.iter().position(|&v| v >= 0.05).unwrap_or(n_pts);

    if start_idx + 2 < n_pts {
        // Second-order forward difference at start_idx
        let dydv_start =
            (-3.0 * ln_j[start_idx] + 4.0 * ln_j[start_idx + 1] - ln_j[start_idx + 2]) / (2.0 * dv);
        if dydv_start > 1e-6 {
            ideality[start_idx] = (1.0 / (vt * dydv_start)).clamp(0.8, 3.0);
        }

        // Central differences for interior forward points
        for i in (start_idx + 1)..(n_pts - 1) {
            let dydv = (ln_j[i + 1] - ln_j[i - 1]) / (2.0 * dv);
            if dydv > 1e-6 {
                ideality[i] = (1.0 / (vt * dydv)).clamp(0.8, 3.0);
            }
        }

        // Second-order backward difference at the final point (V = 1.0 V)
        let last_idx = n_pts - 1;
        let dydv_end =
            (3.0 * ln_j[last_idx] - 4.0 * ln_j[last_idx - 1] + ln_j[last_idx - 2]) / (2.0 * dv);
        if dydv_end > 1e-6 {
            ideality[last_idx] = (1.0 / (vt * dydv_end)).clamp(0.8, 3.0);
        }
    }

    // Operating point evaluation
    let vj_op = solve_vj(v_bias);
    let (jd_op, js_op) = eval_currents(vj_op);
    let j_op = jd_op + js_op;

    // Operating ideality factor using central difference around V_bias
    let delta = 0.005;
    let vj_op_hi = solve_vj(v_bias + delta);
    let (jd_hi, js_hi) = eval_currents(vj_op_hi);
    let j_hi = (jd_hi + js_hi).max(1e-25);

    let vj_op_lo = solve_vj(v_bias - delta);
    let (jd_lo, js_lo) = eval_currents(vj_op_lo);
    let j_lo = (jd_lo + js_lo).max(1e-25);

    let ideality_op = if v_bias >= 0.05 && j_hi > j_lo {
        let d_ln_j = (j_hi / j_lo).ln();
        let dydv = d_ln_j / (2.0 * delta);
        if dydv > 1e-6 {
            (1.0 / (vt * dydv)).clamp(0.8, 3.0)
        } else {
            2.0
        }
    } else {
        2.0
    };

    // Total saturation current J0 (evaluated around V = -0.5 V)
    let (_, js_rev) = eval_currents(-0.5);
    let j0_total = j0_diff + js_rev.abs();

    IvCurveResult {
        v_v,
        j_total_a_cm2,
        j_total_signed_a_cm2,
        j_diff_a_cm2,
        j_srh_a_cm2,
        ideality,
        j0_saturation: j0_total,
        j_operating: j_op,
        ideality_operating: ideality_op,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::materials::MaterialType;

    #[test]
    fn test_monotonic_forward_bias_current() {
        let props = MaterialType::Silicon.properties();
        let iv = compute_iv_curve(&props, &props, 300.0, 1e16, 1e16, 0.5, 1e-7, 0.7);

        // Current must strictly increase for forward bias
        let mut last_j = 0.0;
        for i in 0..iv.v_v.len() {
            if iv.v_v[i] > 0.1 {
                let j = iv.j_total_signed_a_cm2[i];
                assert!(
                    j > last_j,
                    "Current must be monotonically increasing in forward bias"
                );
                last_j = j;
            }
        }
    }

    #[test]
    fn test_ideality_factor_limits() {
        let props = MaterialType::Silicon.properties();
        let iv = compute_iv_curve(&props, &props, 300.0, 1e16, 1e16, 0.3, 1e-7, 0.7);

        // In moderate forward bias (0.1 V to 0.5 V), ideality factor should be in [0.8, 2.2]
        for i in 0..iv.v_v.len() {
            let v = iv.v_v[i];
            if v >= 0.1 && v <= 0.5 {
                let eta = iv.ideality[i];
                assert!(
                    eta >= 0.8 && eta <= 2.2,
                    "Ideality factor at V = {:.2} V was {}, expected in [0.8, 2.2]",
                    v,
                    eta
                );
            }
        }
    }
}
