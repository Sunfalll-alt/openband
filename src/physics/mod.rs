pub mod iv_curve;
pub mod materials;
pub mod metrics;
pub mod solver;

pub use iv_curve::{IvCurveResult, compute_iv_curve};
pub use materials::{EPS_0, K_BOLTZ, MaterialProperties, MaterialType, Q_ELEM, thermal_voltage};
pub use metrics::DeviceMetrics;
pub use solver::{DeviceParams, NUM_NODES, SimulationResult, SolverWorkspace};
