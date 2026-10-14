//! Semiconductor material definitions and temperature-dependent properties.
//!
//! Models temperature-dependent bandgap via the Varshni equation:
//!     Eg(T) = Eg(0) - (alpha * T^2) / (T + beta)
//! and effective densities of states Nc(T), Nv(T) proportional to (T / 300)^1.5.
//! Includes Electron Affinity (chi) for Anderson's heterojunction rule.

use serde::{Deserialize, Serialize};

/// Fundamental physical constants in CGS / semiconductor engineering units.
pub const Q_ELEM: f64 = 1.602_176_634e-19; // Coulombs (C)
pub const K_BOLTZ: f64 = 1.380_649e-23; // Joules / Kelvin (J/K)
pub const EPS_0: f64 = 8.854_187_8128e-14; // Vacuum permittivity in F/cm

/// Computes the thermal voltage Vt = k_B * T / q in Volts.
#[inline]
pub fn thermal_voltage(temp_k: f64) -> f64 {
    (K_BOLTZ * temp_k.max(1.0)) / Q_ELEM
}

/// Supported semiconductor materials.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MaterialType {
    Silicon,
    Germanium,
    GalliumArsenide,
    SiliconCarbide4H,
    AlGaAs,
}

impl MaterialType {
    pub const ALL: [Self; 5] = [
        Self::Silicon,
        Self::Germanium,
        Self::GalliumArsenide,
        Self::SiliconCarbide4H,
        Self::AlGaAs,
    ];

    pub fn display_name(&self) -> &'static str {
        match self {
            Self::Silicon => "Silicon (Si)",
            Self::Germanium => "Germanium (Ge)",
            Self::GalliumArsenide => "Gallium Arsenide (GaAs)",
            Self::SiliconCarbide4H => "Silicon Carbide (4H-SiC)",
            Self::AlGaAs => "AlGaAs (Al0.3Ga0.7As)",
        }
    }

    pub fn short_name(&self) -> &'static str {
        match self {
            Self::Silicon => "Si",
            Self::Germanium => "Ge",
            Self::GalliumArsenide => "GaAs",
            Self::SiliconCarbide4H => "4H-SiC",
            Self::AlGaAs => "Al0.3Ga0.7As",
        }
    }

    pub fn properties(&self) -> MaterialProperties {
        match self {
            Self::Silicon => MaterialProperties {
                material: *self,
                name: "Silicon (Si)",
                eg_0: 1.166,             // eV
                alpha: 4.73e-4,          // eV/K
                beta: 636.0,             // K
                nc_300: 2.80e19,         // cm^-3
                nv_300: 1.04e19,         // cm^-3
                eps_r: 11.7,             // Relative permittivity
                electron_affinity: 4.05, // eV (chi)
                mu_n: 1400.0,            // cm^2/(V*s)
                mu_p: 450.0,             // cm^2/(V*s)
            },
            Self::Germanium => MaterialProperties {
                material: *self,
                name: "Germanium (Ge)",
                eg_0: 0.7437,
                alpha: 4.77e-4,
                beta: 235.0,
                nc_300: 1.04e19,
                nv_300: 6.00e18,
                eps_r: 16.0,
                electron_affinity: 4.00,
                mu_n: 3900.0,
                mu_p: 1900.0,
            },
            Self::GalliumArsenide => MaterialProperties {
                material: *self,
                name: "Gallium Arsenide (GaAs)",
                eg_0: 1.519,
                alpha: 5.405e-4,
                beta: 204.0,
                nc_300: 4.70e17,
                nv_300: 7.00e18,
                eps_r: 12.9,
                electron_affinity: 4.07,
                mu_n: 8500.0,
                mu_p: 400.0,
            },
            Self::SiliconCarbide4H => MaterialProperties {
                material: *self,
                name: "Silicon Carbide (4H-SiC)",
                eg_0: 3.265,
                alpha: 3.30e-4,
                beta: 1200.0,
                nc_300: 1.70e19,
                nv_300: 2.50e19,
                eps_r: 9.7,
                electron_affinity: 3.70, // Standard 3.70 eV
                mu_n: 900.0,
                mu_p: 120.0,
            },
            Self::AlGaAs => MaterialProperties {
                material: *self,
                name: "AlGaAs (Al0.3Ga0.7As)",
                eg_0: 1.8945, // Eg(300K) ~ 1.798 eV
                alpha: 5.405e-4,
                beta: 204.0,
                nc_300: 6.50e17, // cm^-3
                nv_300: 1.10e19, // cm^-3
                eps_r: 12.05,
                electron_affinity: 3.74, // eV
                mu_n: 2500.0,
                mu_p: 150.0,
            },
        }
    }
}

/// Physical parameters of a semiconductor material.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct MaterialProperties {
    pub material: MaterialType,
    pub name: &'static str,
    /// Bandgap at 0 Kelvin (eV)
    pub eg_0: f64,
    /// Varshni parameter alpha (eV/K)
    pub alpha: f64,
    /// Varshni parameter beta (K)
    pub beta: f64,
    /// Effective conduction band DOS at 300 K (cm^-3)
    pub nc_300: f64,
    /// Effective valence band DOS at 300 K (cm^-3)
    pub nv_300: f64,
    /// Relative dielectric permittivity epsilon_r
    pub eps_r: f64,
    /// Electron affinity chi (eV)
    pub electron_affinity: f64,
    /// Electron low-field mobility (cm^2 / (V*s))
    pub mu_n: f64,
    /// Hole low-field mobility (cm^2 / (V*s))
    pub mu_p: f64,
}

impl MaterialProperties {
    /// Computes bandgap Eg(T) in eV using the Varshni empirical relation.
    #[inline]
    pub fn bandgap(&self, temp_k: f64) -> f64 {
        let t = temp_k.max(1.0);
        let eg = self.eg_0 - (self.alpha * t * t) / (t + self.beta);
        eg.max(0.05)
    }

    /// Effective density of states in the conduction band Nc(T) in cm^-3.
    #[inline]
    pub fn nc(&self, temp_k: f64) -> f64 {
        let t = temp_k.max(1.0);
        self.nc_300 * (t / 300.0).powf(1.5)
    }

    /// Effective density of states in the valence band Nv(T) in cm^-3.
    #[inline]
    pub fn nv(&self, temp_k: f64) -> f64 {
        let t = temp_k.max(1.0);
        self.nv_300 * (t / 300.0).powf(1.5)
    }

    /// Intrinsic carrier concentration n_i(T) in cm^-3.
    /// Clamped to avoid floating point underflow for wide bandgaps at cryogenic temperatures.
    #[inline]
    pub fn ni(&self, temp_k: f64) -> f64 {
        let t = temp_k.max(1.0);
        let vt = thermal_voltage(t);
        let eg = self.bandgap(t);
        let nc = self.nc(t);
        let nv = self.nv(t);
        let exponent = -eg / (2.0 * vt);

        if exponent < -700.0 {
            1e-35
        } else {
            let ni_val = (nc * nv).sqrt() * exponent.exp();
            ni_val.max(1e-35)
        }
    }

    /// Intrinsic Debye length L_D = sqrt(eps_0 * eps_r * Vt / (q * n_i)) in cm.
    #[inline]
    pub fn debye_length(&self, temp_k: f64) -> f64 {
        let vt = thermal_voltage(temp_k);
        let ni = self.ni(temp_k);
        let numerator = EPS_0 * self.eps_r * vt;
        let denominator = Q_ELEM * ni;
        (numerator / denominator).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_silicon_at_300k() {
        let props = MaterialType::Silicon.properties();
        let eg = props.bandgap(300.0);
        assert!(
            (eg - 1.12).abs() < 0.02,
            "Silicon Eg at 300K should be ~1.12 eV, got {}",
            eg
        );

        let ni = props.ni(300.0);
        assert!(
            ni > 1e9 && ni < 5e10,
            "Silicon ni at 300K should be around 1e10, got {}",
            ni
        );

        let vt = thermal_voltage(300.0);
        assert!(
            (vt - 0.02585).abs() < 0.001,
            "Thermal voltage at 300K should be ~25.85 mV, got {}",
            vt
        );
        assert_eq!(props.electron_affinity, 4.05);
    }

    #[test]
    fn test_algaas_at_300k() {
        let props = MaterialType::AlGaAs.properties();
        let eg = props.bandgap(300.0);
        assert!(
            (eg - 1.798).abs() < 0.01,
            "AlGaAs Eg at 300K should be ~1.798 eV, got {}",
            eg
        );
        assert_eq!(props.eps_r, 12.05);
        assert_eq!(props.electron_affinity, 3.74);
    }

    #[test]
    fn test_all_materials_positive_bandgap() {
        for mat in MaterialType::ALL {
            let props = mat.properties();
            for t in [100.0, 300.0, 600.0] {
                let eg = props.bandgap(t);
                assert!(
                    eg > 0.0,
                    "Bandgap for {:?} at {} K must be positive",
                    mat,
                    t
                );
                let ni = props.ni(t);
                assert!(ni > 0.0, "ni for {:?} at {} K must be positive", mat, t);
            }
        }
    }
}
