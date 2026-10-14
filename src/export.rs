//! CSV export functionality for simulation profiles and device metadata.

use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::physics::{DeviceParams, SimulationResult};

/// Generates a standardized, comma-separated values (CSV) string from a simulation result.
pub fn generate_csv_string(params: &DeviceParams, result: &SimulationResult) -> String {
    let mut out = String::with_capacity(32 * 1024);

    // Metadata header
    out.push_str("# openband 1D Semiconductor Simulation Profile\n");
    if params.is_heterojunction {
        out.push_str(&format!(
            "# Junction Type: Heterojunction (p-side: {} | n-side: {})\n",
            params.material.display_name(),
            params.material_n.display_name()
        ));
        out.push_str(&format!(
            "# Band Offsets: Delta Ec = {:.3} eV | Delta Ev = {:.3} eV\n",
            result.metrics.delta_ec_ev, result.metrics.delta_ev_ev
        ));
    } else {
        out.push_str(&format!(
            "# Junction Type: Homojunction ({})\n",
            params.material.display_name()
        ));
    }
    out.push_str(&format!("# Temperature: {:.2} K\n", params.temp_k));
    out.push_str(&format!(
        "# Carrier Lifetime tau_0: {:.3e} s\n",
        params.tau_0
    ));
    out.push_str(&format!("# Acceptor Doping Na: {:.3e} cm^-3\n", params.na));
    out.push_str(&format!("# Donor Doping Nd: {:.3e} cm^-3\n", params.nd));
    out.push_str(&format!("# Device Length: {:.3} um\n", params.length_um));
    out.push_str(&format!(
        "# Junction Position: {:.3} um\n",
        params.junction_ratio * params.length_um
    ));
    out.push_str(&format!("# Applied Bias: {:.3} V\n", params.v_bias));
    out.push_str(&format!(
        "# Analytical V_bi: {:.4} V | Numerical V_bi: {:.4} V\n",
        result.metrics.v_bi_analytical, result.metrics.v_bi_numerical
    ));
    out.push_str(&format!(
        "# Analytical W_dep: {:.4} um | Numerical W_dep: {:.4} um\n",
        result.metrics.w_dep_analytical_um, result.metrics.w_dep_numerical_um
    ));
    out.push_str(&format!(
        "# Analytical E_max: {:.2} kV/cm | Numerical E_max: {:.2} kV/cm\n",
        result.metrics.e_max_analytical_kv_cm, result.metrics.e_max_numerical_kv_cm
    ));
    out.push_str(&format!(
        "# Analytical C_j: {:.4} nF/cm^2 | Numerical C_j: {:.4} nF/cm^2\n",
        result.metrics.c_j_analytical_nf_cm2, result.metrics.c_j_numerical_nf_cm2
    ));
    out.push_str(&format!(
        "# Diode J0 Saturation: {:.3e} A/cm^2 | J(V_bias): {:.3e} A/cm^2 | Ideality eta: {:.2}\n",
        result.iv_curve.j0_saturation,
        result.iv_curve.j_operating,
        result.iv_curve.ideality_operating
    ));
    out.push_str("#\n");

    // Table Column Headers
    out.push_str(
        "x_um,psi_V,e_field_kV_cm,rho_C_cm3,n_cm3,p_cm3,Ec_eV,Ev_eV,Ei_eV,EFn_eV,EFp_eV,net_doping_cm3\n",
    );

    let n = result.x_um.len();
    for i in 0..n {
        out.push_str(&format!(
            "{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}\n",
            result.x_um[i],
            result.psi_v[i],
            result.e_field_kv_cm[i],
            result.rho_c_cm3[i],
            result.n_cm3[i],
            result.p_cm3[i],
            result.e_c_ev[i],
            result.e_v_ev[i],
            result.e_i_ev[i],
            result.e_fn_ev[i],
            result.e_fp_ev[i],
            result.net_doping_cm3[i],
        ));
    }

    out
}

/// Writes the simulation profile CSV to a designated filesystem path.
pub fn export_to_file<P: AsRef<Path>>(
    params: &DeviceParams,
    result: &SimulationResult,
    path: P,
) -> std::io::Result<()> {
    let content = generate_csv_string(params, result);
    let mut file = File::create(path)?;
    file.write_all(content.as_bytes())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::solver::solve;

    #[test]
    fn test_csv_generation_format() {
        let params = DeviceParams::default();
        let res = solve(&params);
        let csv = generate_csv_string(&params, &res);

        assert!(csv.contains("# openband 1D Semiconductor Simulation Profile"));
        assert!(csv.contains("x_um,psi_V,e_field_kV_cm"));
        let lines: Vec<&str> = csv.lines().collect();
        // 15 metadata lines + header line + 600 data lines = 616 lines
        assert_eq!(lines.len(), 616);
    }
}
