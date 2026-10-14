//! Main GUI Application: State, Parameter Inspector, and Scientific Plot Viewports.

use std::time::Instant;

use eframe::egui::{self, Color32, RichText, ScrollArea, Vec2};
use egui_plot::{Corner, HLine, Legend, Line, LineStyle, Plot, PlotPoints, VLine};

use crate::export::export_to_file;
use crate::physics::{
    DeviceParams, NUM_NODES, SimulationResult, SolverWorkspace, materials::MaterialType,
};

/// Active plot viewport arrangement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlotViewMode {
    QuadSplit,
    EnergyBands,
    CarrierConcentrations,
    ElectricFieldAndPotential,
    SpaceChargeDensity,
    IvCharacteristic,
}

impl PlotViewMode {
    pub const ALL: [Self; 6] = [
        Self::QuadSplit,
        Self::EnergyBands,
        Self::CarrierConcentrations,
        Self::ElectricFieldAndPotential,
        Self::SpaceChargeDensity,
        Self::IvCharacteristic,
    ];

    pub fn display_name(&self) -> &'static str {
        match self {
            Self::QuadSplit => "⊞ Quad Split View (2x2 Grid)",
            Self::EnergyBands => "⚡ Energy Bands (eV)",
            Self::CarrierConcentrations => "📊 Carrier Concentrations (log10 cm^-3)",
            Self::ElectricFieldAndPotential => "⚡ Field (kV/cm) & Potential (V)",
            Self::SpaceChargeDensity => "🔋 Space Charge Density (C/cm^3)",
            Self::IvCharacteristic => "📈 I-V Characteristic J(V)",
        }
    }
}

/// Display mode for Electric Field and Potential scaling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldPotentialMode {
    OverlaidScaled, // E(x) [kV/cm] + Potential psi(x) scaled x10 (0.1 V/div)
    Separated,      // Two separate stacked subplots
}

/// The main openband application state.
pub struct OpenBandApp {
    pub params: DeviceParams,
    pub log_na: f64,
    pub log_nd: f64,
    pub log_tau_0: f64,
    pub view_mode: PlotViewMode,
    pub field_potential_mode: FieldPotentialMode,
    pub show_legends: bool,
    pub simulation: SimulationResult,
    pub workspace: SolverWorkspace,
    pub last_solved_params: Option<DeviceParams>,
    pub status_message: Option<(String, Instant)>,
}

impl Default for OpenBandApp {
    fn default() -> Self {
        let params = DeviceParams::default();
        let log_na = params.na.log10();
        let log_nd = params.nd.log10();
        let log_tau_0 = params.tau_0.log10();
        let mut workspace = SolverWorkspace::new(NUM_NODES);
        let simulation = workspace.solve(&params);
        let last_solved_params = Some(params.clone());

        Self {
            params,
            log_na,
            log_nd,
            log_tau_0,
            view_mode: PlotViewMode::QuadSplit,
            field_potential_mode: FieldPotentialMode::OverlaidScaled,
            show_legends: true,
            simulation,
            workspace,
            last_solved_params,
            status_message: None,
        }
    }
}

impl OpenBandApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self::default()
    }

    /// Re-evaluates Poisson-Boltzmann if parameters were modified.
    pub fn update_simulation_if_needed(&mut self) {
        let mut need_solve = false;
        if let Some(ref last) = self.last_solved_params {
            if last.material != self.params.material
                || last.material_n != self.params.material_n
                || last.is_heterojunction != self.params.is_heterojunction
                || (last.temp_k - self.params.temp_k).abs() > 1e-6
                || (last.na - self.params.na).abs() > 1e-6
                || (last.nd - self.params.nd).abs() > 1e-6
                || (last.length_um - self.params.length_um).abs() > 1e-6
                || (last.junction_ratio - self.params.junction_ratio).abs() > 1e-6
                || (last.v_bias - self.params.v_bias).abs() > 1e-6
                || (last.tau_0 - self.params.tau_0).abs() > 1e-12
            {
                need_solve = true;
            }
        } else {
            need_solve = true;
        }

        if need_solve {
            self.simulation = self.workspace.solve(&self.params);
            self.last_solved_params = Some(self.params.clone());
        }
    }

    /// Sets device parameters and syncs log sliders.
    pub fn apply_preset(&mut self, preset: DeviceParams) {
        self.params = preset;
        self.log_na = self.params.na.log10();
        self.log_nd = self.params.nd.log10();
        self.log_tau_0 = self.params.tau_0.log10();
        self.update_simulation_if_needed();
    }
}

impl eframe::App for OpenBandApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Automatically check if any parameter changed and re-solve in < 0.2 ms
        self.update_simulation_if_needed();

        // 1. Left SidePanel: Parameter Inspector & Telemetry
        egui::SidePanel::left("parameter_inspector_panel")
            .resizable(true)
            .default_width(360.0)
            .min_width(320.0)
            .max_width(480.0)
            .show(ctx, |ui| {
                ScrollArea::vertical().show(ui, |ui| {
                    ui.add_space(8.0);
                    ui.heading(RichText::new("openband v0.2.0").strong());
                    ui.label(
                        RichText::new("1D Semiconductor & Heterojunction Simulator")
                            .small()
                            .weak(),
                    );
                    ui.separator();

                    // --- Heterojunction Mode Toggle ---
                    ui.checkbox(
                        &mut self.params.is_heterojunction,
                        RichText::new("Heterojunction Mode (Separate p/n Materials)").strong(),
                    );
                    ui.add_space(4.0);

                    if self.params.is_heterojunction {
                        // Heterojunction Mode: separate p-side and n-side dropdowns
                        ui.label(RichText::new("p-side Material (x < xj):").strong());
                        egui::ComboBox::from_id_salt("material_p_combo")
                            .selected_text(self.params.material.display_name())
                            .show_ui(ui, |ui| {
                                for mat in MaterialType::ALL {
                                    ui.selectable_value(
                                        &mut self.params.material,
                                        mat,
                                        mat.display_name(),
                                    );
                                }
                            });
                        let props_p = self.params.material.properties();
                        let eg_p = props_p.bandgap(self.params.temp_k);
                        ui.label(
                            RichText::new(format!(
                                "p-side: Eg: {:.3} eV | chi: {:.2} eV | eps_r: {:.1}",
                                eg_p, props_p.electron_affinity, props_p.eps_r
                            ))
                            .small()
                            .weak(),
                        );

                        ui.add_space(4.0);
                        ui.label(RichText::new("n-side Material (x >= xj):").strong());
                        egui::ComboBox::from_id_salt("material_n_combo")
                            .selected_text(self.params.material_n.display_name())
                            .show_ui(ui, |ui| {
                                for mat in MaterialType::ALL {
                                    ui.selectable_value(
                                        &mut self.params.material_n,
                                        mat,
                                        mat.display_name(),
                                    );
                                }
                            });
                        let props_n = self.params.material_n.properties();
                        let eg_n = props_n.bandgap(self.params.temp_k);
                        ui.label(
                            RichText::new(format!(
                                "n-side: Eg: {:.3} eV | chi: {:.2} eV | eps_r: {:.1}",
                                eg_n, props_n.electron_affinity, props_n.eps_r
                            ))
                            .small()
                            .weak(),
                        );
                    } else {
                        // Homojunction Mode
                        ui.label(RichText::new("Semiconductor Material").strong());
                        egui::ComboBox::from_id_salt("material_combo")
                            .selected_text(self.params.material.display_name())
                            .show_ui(ui, |ui| {
                                for mat in MaterialType::ALL {
                                    ui.selectable_value(
                                        &mut self.params.material,
                                        mat,
                                        mat.display_name(),
                                    );
                                }
                            });
                        let props = self.params.material.properties();
                        let eg = props.bandgap(self.params.temp_k);
                        let ni = props.ni(self.params.temp_k);
                        ui.label(
                            RichText::new(format!(
                                "Eg(T): {:.3} eV | chi: {:.2} eV | eps_r: {:.1} | ni: {:.2e} cm^-3",
                                eg, props.electron_affinity, props.eps_r, ni
                            ))
                            .small()
                            .weak(),
                        );
                    }
                    ui.add_space(6.0);
                    ui.separator();

                    // --- Doping Concentrations ---
                    ui.label(RichText::new("Doping Concentrations (Logarithmic)").strong());
                    let prev_log_na = self.log_na;
                    ui.add(
                        egui::Slider::new(&mut self.log_na, 14.0..=19.0)
                            .text("log10(Na) [cm^-3]")
                            .step_by(0.01),
                    );
                    if (self.log_na - prev_log_na).abs() > 1e-9 {
                        self.params.na = 10f64.powf(self.log_na);
                    }
                    ui.label(
                        RichText::new(format!("Na (p-side): {:.2e} cm^-3", self.params.na))
                            .small()
                            .color(Color32::from_rgb(255, 170, 90)),
                    );

                    ui.add_space(4.0);
                    let prev_log_nd = self.log_nd;
                    ui.add(
                        egui::Slider::new(&mut self.log_nd, 14.0..=19.0)
                            .text("log10(Nd) [cm^-3]")
                            .step_by(0.01),
                    );
                    if (self.log_nd - prev_log_nd).abs() > 1e-9 {
                        self.params.nd = 10f64.powf(self.log_nd);
                    }
                    ui.label(
                        RichText::new(format!("Nd (n-side): {:.2e} cm^-3", self.params.nd))
                            .small()
                            .color(Color32::from_rgb(80, 210, 255)),
                    );
                    ui.add_space(6.0);
                    ui.separator();

                    // --- Operating & Geometry Controls ---
                    ui.label(RichText::new("Operating & Geometry Parameters").strong());

                    ui.add(
                        egui::Slider::new(&mut self.params.temp_k, 100.0..=600.0)
                            .text("Temperature T (K)")
                            .step_by(1.0),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.v_bias, -5.0..=0.8)
                            .text("Applied Bias V_bias (V)")
                            .step_by(0.01),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.length_um, 0.2..=5.0)
                            .text("Device Length L (um)")
                            .step_by(0.05),
                    );

                    ui.add(
                        egui::Slider::new(&mut self.params.junction_ratio, 0.1..=0.9)
                            .text("Junction xj / L")
                            .step_by(0.01),
                    );
                    let xj_um = self.params.junction_ratio * self.params.length_um;
                    ui.label(
                        RichText::new(format!("Junction Position xj: {:.3} um", xj_um))
                            .small()
                            .weak(),
                    );

                    ui.add_space(4.0);
                    ui.label(RichText::new("Recombination Lifetime tau_0:").strong());
                    let prev_log_tau = self.log_tau_0;
                    ui.add(
                        egui::Slider::new(&mut self.log_tau_0, -9.0..=-5.0)
                            .text("log10(tau_0) [s]")
                            .step_by(0.05),
                    );
                    if (self.log_tau_0 - prev_log_tau).abs() > 1e-9 {
                        self.params.tau_0 = 10f64.powf(self.log_tau_0);
                    }
                    ui.label(
                        RichText::new(format!(
                            "tau_0: {:.1e} s ({:.1} ns)",
                            self.params.tau_0,
                            self.params.tau_0 * 1e9
                        ))
                        .small()
                        .weak(),
                    );

                    ui.add_space(6.0);
                    ui.separator();

                    // --- Quick Presets ---
                    ui.label(RichText::new("Quick Physics Presets").strong());
                    ui.horizontal_wrapped(|ui| {
                        if ui.button("Symmetric Si Junction").clicked() {
                            self.apply_preset(DeviceParams {
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
                            });
                        }
                        if ui.button("Asymmetric p+n Diode").clicked() {
                            self.apply_preset(DeviceParams {
                                material: MaterialType::Silicon,
                                material_n: MaterialType::Silicon,
                                is_heterojunction: false,
                                temp_k: 300.0,
                                na: 1.0e18,
                                nd: 1.0e15,
                                length_um: 2.0,
                                junction_ratio: 0.3,
                                v_bias: 0.0,
                                tau_0: 1.0e-7,
                            });
                        }
                        if ui.button("4H-SiC Power Diode").clicked() {
                            self.apply_preset(DeviceParams {
                                material: MaterialType::SiliconCarbide4H,
                                material_n: MaterialType::SiliconCarbide4H,
                                is_heterojunction: false,
                                temp_k: 300.0,
                                na: 1.0e18,
                                nd: 1.0e15,
                                length_um: 4.0,
                                junction_ratio: 0.2,
                                v_bias: 0.0,
                                tau_0: 1.0e-7,
                            });
                        }
                        if ui.button("AlGaAs / GaAs Heterojunction").clicked() {
                            self.apply_preset(DeviceParams {
                                material: MaterialType::AlGaAs,
                                material_n: MaterialType::GalliumArsenide,
                                is_heterojunction: true,
                                temp_k: 300.0,
                                na: 5.0e17,
                                nd: 1.0e16,
                                length_um: 1.0,
                                junction_ratio: 0.4,
                                v_bias: 0.0,
                                tau_0: 1.0e-7,
                            });
                        }
                        if ui.button("Cryogenic 100K Mode").clicked() {
                            self.apply_preset(DeviceParams {
                                material: MaterialType::Silicon,
                                material_n: MaterialType::Silicon,
                                is_heterojunction: false,
                                temp_k: 100.0,
                                na: 1.0e16,
                                nd: 1.0e16,
                                length_um: 1.0,
                                junction_ratio: 0.5,
                                v_bias: 0.0,
                                tau_0: 1.0e-7,
                            });
                        }
                    });

                    ui.add_space(8.0);
                    if ui
                        .button(RichText::new("💾 Export Simulation to CSV").strong())
                        .clicked()
                    {
                        let path = "openband_simulation.csv";
                        match export_to_file(&self.params, &self.simulation, path) {
                            Ok(()) => {
                                self.status_message = Some((
                                    format!("Exported 600 nodes to {}", path),
                                    Instant::now(),
                                ));
                            }
                            Err(e) => {
                                self.status_message =
                                    Some((format!("Export failed: {}", e), Instant::now()));
                            }
                        }
                    }

                    if let Some((ref msg, time)) = self.status_message {
                        if time.elapsed().as_secs() < 5 {
                            ui.label(RichText::new(msg).color(Color32::GREEN).small());
                        }
                    }

                    ui.add_space(8.0);
                    ui.separator();

                    // --- Live Telemetry Readout Table ---
                    ui.label(RichText::new("Analytical vs. Numerical Telemetry").strong());
                    let m = &self.simulation.metrics;

                    egui::Grid::new("telemetry_table")
                        .striped(true)
                        .min_col_width(80.0)
                        .show(ui, |ui| {
                            ui.label(RichText::new("Parameter").strong());
                            ui.label(RichText::new("Depletion (Ana.)").strong());
                            ui.label(RichText::new("Poisson-Boltzmann").strong());
                            ui.end_row();

                            ui.label("V_bi");
                            ui.label(format!("{:.4} V", m.v_bi_analytical));
                            ui.label(format!("{:.4} V", m.v_bi_numerical));
                            ui.end_row();

                            ui.label("W_dep");
                            ui.label(format!("{:.4} um", m.w_dep_analytical_um));
                            ui.label(format!("{:.4} um", m.w_dep_numerical_um));
                            ui.end_row();

                            ui.label("|E_max|");
                            ui.label(format!("{:.2} kV/cm", m.e_max_analytical_kv_cm));
                            ui.label(format!("{:.2} kV/cm", m.e_max_numerical_kv_cm));
                            ui.end_row();

                            ui.label("C_j");
                            ui.label(format!("{:.3} nF/cm^2", m.c_j_analytical_nf_cm2));
                            ui.label(format!("{:.3} nF/cm^2", m.c_j_numerical_nf_cm2));
                            ui.end_row();
                        });

                    if self.params.is_heterojunction {
                        ui.add_space(6.0);
                        ui.label(RichText::new("Heterojunction Band Offsets (Anderson)").strong());
                        ui.label(format!("Delta Ec (Conduction): {:.3} eV", m.delta_ec_ev));
                        ui.label(format!("Delta Ev (Valence): {:.3} eV", m.delta_ev_ev));
                    }

                    ui.add_space(6.0);
                    ui.label(RichText::new("Diode Transport & I-V Telemetry").strong());
                    let iv = &self.simulation.iv_curve;
                    ui.label(format!("Saturation J0: {:.3e} A/cm^2", iv.j0_saturation));
                    ui.label(format!("J(V_bias): {:.3e} A/cm^2", iv.j_operating));
                    ui.label(format!("Ideality Factor eta: {:.2}", iv.ideality_operating));

                    ui.add_space(8.0);
                    ui.group(|ui| {
                        let status_text = if self.simulation.converged {
                            RichText::new("● FDM Solver Converged").color(Color32::GREEN)
                        } else {
                            RichText::new("▲ Iteration limit reached").color(Color32::YELLOW)
                        };
                        ui.label(status_text);
                        ui.label(
                            RichText::new(format!(
                                "Iterations: {} | Final max|δu|: {:.2e}",
                                self.simulation.iterations, self.simulation.max_delta_u
                            ))
                            .small()
                            .weak(),
                        );
                    });

                    ui.add_space(12.0);
                });
            });

        // 2. Top Panel: View Mode Navigation & Toolbar Controls
        egui::TopBottomPanel::top("top_navigation_bar").show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("Plot Arrangement:").strong());
                for mode in PlotViewMode::ALL {
                    let is_active = self.view_mode == mode;
                    if ui
                        .selectable_label(is_active, mode.display_name())
                        .clicked()
                    {
                        self.view_mode = mode;
                    }
                }

                ui.separator();
                ui.checkbox(&mut self.show_legends, "Show Plot Legends");

                if self.view_mode == PlotViewMode::ElectricFieldAndPotential
                    || self.view_mode == PlotViewMode::QuadSplit
                {
                    ui.separator();
                    ui.label(RichText::new("Field & Potential:").small());
                    let is_overlaid =
                        self.field_potential_mode == FieldPotentialMode::OverlaidScaled;
                    if ui
                        .selectable_label(is_overlaid, "Overlaid (x10 V)")
                        .clicked()
                    {
                        self.field_potential_mode = FieldPotentialMode::OverlaidScaled;
                    }
                    let is_separated = self.field_potential_mode == FieldPotentialMode::Separated;
                    if ui
                        .selectable_label(is_separated, "Separated Subplots")
                        .clicked()
                    {
                        self.field_potential_mode = FieldPotentialMode::Separated;
                    }
                }
            });
        });

        // 3. CentralPanel: Interactive Scientific Plots
        egui::CentralPanel::default().show(ctx, |ui| {
            let xj_um = self.params.junction_ratio * self.params.length_um;
            let res = &self.simulation;
            let show_legends = self.show_legends;
            let fp_mode = self.field_potential_mode;

            match self.view_mode {
                PlotViewMode::QuadSplit => {
                    let total_size = ui.available_size();
                    let half_w = (total_size.x - 16.0) * 0.5;
                    let half_h = (total_size.y - 16.0) * 0.5;

                    ui.vertical(|ui| {
                        ui.horizontal(|ui| {
                            ui.allocate_ui(Vec2::new(half_w, half_h), |ui| {
                                Self::draw_energy_bands(ui, res, xj_um, show_legends);
                            });
                            ui.allocate_ui(Vec2::new(half_w, half_h), |ui| {
                                Self::draw_carrier_concentrations(ui, res, xj_um, show_legends);
                            });
                        });
                        ui.add_space(8.0);
                        ui.horizontal(|ui| {
                            ui.allocate_ui(Vec2::new(half_w, half_h), |ui| {
                                Self::draw_electric_field_and_potential(
                                    ui,
                                    res,
                                    xj_um,
                                    show_legends,
                                    fp_mode,
                                );
                            });
                            ui.allocate_ui(Vec2::new(half_w, half_h), |ui| {
                                Self::draw_space_charge_density(ui, res, xj_um, show_legends);
                            });
                        });
                    });
                }
                PlotViewMode::EnergyBands => {
                    Self::draw_energy_bands(ui, res, xj_um, show_legends);
                }
                PlotViewMode::CarrierConcentrations => {
                    Self::draw_carrier_concentrations(ui, res, xj_um, show_legends);
                }
                PlotViewMode::ElectricFieldAndPotential => {
                    Self::draw_electric_field_and_potential(ui, res, xj_um, show_legends, fp_mode);
                }
                PlotViewMode::SpaceChargeDensity => {
                    Self::draw_space_charge_density(ui, res, xj_um, show_legends);
                }
                PlotViewMode::IvCharacteristic => {
                    Self::draw_iv_characteristic(ui, res, self.params.v_bias, show_legends);
                }
            }
        });
    }
}

impl OpenBandApp {
    /// Renders Energy Bands: Ec, Ev, Ei, EFn, EFp
    fn draw_energy_bands(
        ui: &mut egui::Ui,
        res: &SimulationResult,
        xj_um: f64,
        show_legends: bool,
    ) {
        let n = res.x_um.len();
        let pts_ec: Vec<[f64; 2]> = (0..n).map(|i| [res.x_um[i], res.e_c_ev[i]]).collect();
        let pts_ev: Vec<[f64; 2]> = (0..n).map(|i| [res.x_um[i], res.e_v_ev[i]]).collect();
        let pts_ei: Vec<[f64; 2]> = (0..n).map(|i| [res.x_um[i], res.e_i_ev[i]]).collect();
        let pts_efn: Vec<[f64; 2]> = (0..n).map(|i| [res.x_um[i], res.e_fn_ev[i]]).collect();
        let pts_efp: Vec<[f64; 2]> = (0..n).map(|i| [res.x_um[i], res.e_fp_ev[i]]).collect();

        let mut plot = Plot::new("plot_energy_bands")
            .x_axis_label("Position x (um)")
            .y_axis_label("Energy (eV)");

        if show_legends {
            plot = plot.legend(Legend::default().position(Corner::LeftBottom));
        }

        plot.show(ui, |plot_ui| {
            plot_ui.line(
                Line::new(PlotPoints::new(pts_ec))
                    .name("Ec (Conduction Band)")
                    .color(Color32::from_rgb(0, 180, 255))
                    .width(2.0_f32),
            );
            plot_ui.line(
                Line::new(PlotPoints::new(pts_ev))
                    .name("Ev (Valence Band)")
                    .color(Color32::from_rgb(255, 90, 80))
                    .width(2.0_f32),
            );
            plot_ui.line(
                Line::new(PlotPoints::new(pts_ei))
                    .name("Ei (Intrinsic Fermi)")
                    .color(Color32::from_rgb(100, 230, 120))
                    .style(LineStyle::Dashed { length: 5.0 })
                    .width(1.5_f32),
            );
            plot_ui.line(
                Line::new(PlotPoints::new(pts_efn))
                    .name("EFn (Electron Quasi-Fermi)")
                    .color(Color32::from_rgb(255, 200, 50))
                    .style(LineStyle::Dashed { length: 8.0 })
                    .width(1.5_f32),
            );
            plot_ui.line(
                Line::new(PlotPoints::new(pts_efp))
                    .name("EFp (Hole Quasi-Fermi)")
                    .color(Color32::from_rgb(220, 130, 255))
                    .style(LineStyle::Dashed { length: 8.0 })
                    .width(1.5_f32),
            );
            plot_ui.vline(
                VLine::new(xj_um)
                    .name("Metallurgical Junction xj")
                    .color(Color32::from_gray(200))
                    .style(LineStyle::Dashed { length: 6.0 }),
            );
        });
    }

    /// Renders Carrier Concentrations: log10(n), log10(p), log10(Na), log10(Nd)
    fn draw_carrier_concentrations(
        ui: &mut egui::Ui,
        res: &SimulationResult,
        xj_um: f64,
        show_legends: bool,
    ) {
        let n = res.x_um.len();
        let pts_n: Vec<[f64; 2]> = (0..n)
            .map(|i| [res.x_um[i], res.n_cm3[i].max(1e-12).log10()])
            .collect();
        let pts_p: Vec<[f64; 2]> = (0..n)
            .map(|i| [res.x_um[i], res.p_cm3[i].max(1e-12).log10()])
            .collect();
        let pts_doping: Vec<[f64; 2]> = (0..n)
            .map(|i| [res.x_um[i], res.net_doping_cm3[i].abs().max(1e-12).log10()])
            .collect();

        let mut plot = Plot::new("plot_carrier_concentrations")
            .x_axis_label("Position x (um)")
            .y_axis_label("log10(Concentration [cm^-3])");

        if show_legends {
            plot = plot.legend(Legend::default().position(Corner::LeftBottom));
        }

        plot.show(ui, |plot_ui| {
            plot_ui.line(
                Line::new(PlotPoints::new(pts_n))
                    .name("Electrons log10(n)")
                    .color(Color32::from_rgb(60, 160, 255))
                    .width(2.0_f32),
            );
            plot_ui.line(
                Line::new(PlotPoints::new(pts_p))
                    .name("Holes log10(p)")
                    .color(Color32::from_rgb(255, 70, 70))
                    .width(2.0_f32),
            );
            plot_ui.line(
                Line::new(PlotPoints::new(pts_doping))
                    .name("Net Doping log10|Nd - Na|")
                    .color(Color32::from_rgb(200, 200, 200))
                    .style(LineStyle::Dashed { length: 4.0 })
                    .width(1.5_f32),
            );
            plot_ui.vline(
                VLine::new(xj_um)
                    .name("Metallurgical Junction xj")
                    .color(Color32::from_gray(200))
                    .style(LineStyle::Dashed { length: 6.0 }),
            );
        });
    }

    /// Renders Electric Field E(x) [kV/cm] and Potential psi(x) [V] with dual-scale or separated mode
    fn draw_electric_field_and_potential(
        ui: &mut egui::Ui,
        res: &SimulationResult,
        xj_um: f64,
        show_legends: bool,
        mode: FieldPotentialMode,
    ) {
        let n = res.x_um.len();
        let pts_e: Vec<[f64; 2]> = (0..n)
            .map(|i| [res.x_um[i], res.e_field_kv_cm[i]])
            .collect();
        let pts_psi: Vec<[f64; 2]> = (0..n).map(|i| [res.x_um[i], res.psi_v[i]]).collect();

        match mode {
            FieldPotentialMode::OverlaidScaled => {
                // Scale psi(x) by 10 (0.1 V/div) so it is clearly visible alongside kV/cm
                let pts_psi_scaled: Vec<[f64; 2]> =
                    (0..n).map(|i| [res.x_um[i], res.psi_v[i] * 10.0]).collect();

                let mut plot = Plot::new("plot_electric_field_potential_overlaid")
                    .x_axis_label("Position x (um)")
                    .y_axis_label("E (kV/cm) / Potential psi (x10 V)");

                if show_legends {
                    plot = plot.legend(Legend::default().position(Corner::LeftBottom));
                }

                plot.show(ui, |plot_ui| {
                    plot_ui.line(
                        Line::new(PlotPoints::new(pts_e))
                            .name("Electric Field E (kV/cm)")
                            .color(Color32::from_rgb(255, 205, 40))
                            .width(2.0_f32),
                    );
                    plot_ui.line(
                        Line::new(PlotPoints::new(pts_psi_scaled))
                            .name("Potential psi (x10 V)")
                            .color(Color32::from_rgb(175, 115, 255))
                            .width(2.0_f32),
                    );
                    plot_ui.hline(
                        HLine::new(0.0)
                            .color(Color32::from_gray(100))
                            .style(LineStyle::Dashed { length: 5.0 }),
                    );
                    plot_ui.vline(
                        VLine::new(xj_um)
                            .name("Metallurgical Junction xj")
                            .color(Color32::from_gray(200))
                            .style(LineStyle::Dashed { length: 6.0 }),
                    );
                });
            }
            FieldPotentialMode::Separated => {
                let avail_h = ui.available_height();
                let sub_h = (avail_h - 10.0) * 0.5;

                ui.vertical(|ui| {
                    ui.allocate_ui(Vec2::new(ui.available_width(), sub_h), |ui| {
                        let mut plot_e = Plot::new("plot_separated_e")
                            .x_axis_label("Position x (um)")
                            .y_axis_label("Electric Field E (kV/cm)");
                        if show_legends {
                            plot_e = plot_e.legend(Legend::default().position(Corner::LeftBottom));
                        }
                        plot_e.show(ui, |plot_ui| {
                            plot_ui.line(
                                Line::new(PlotPoints::new(pts_e))
                                    .name("Electric Field E (kV/cm)")
                                    .color(Color32::from_rgb(255, 205, 40))
                                    .width(2.0_f32),
                            );
                            plot_ui.vline(
                                VLine::new(xj_um)
                                    .name("Metallurgical Junction xj")
                                    .color(Color32::from_gray(200))
                                    .style(LineStyle::Dashed { length: 6.0 }),
                            );
                        });
                    });

                    ui.add_space(8.0);
                    ui.allocate_ui(Vec2::new(ui.available_width(), sub_h), |ui| {
                        let mut plot_psi = Plot::new("plot_separated_psi")
                            .x_axis_label("Position x (um)")
                            .y_axis_label("Potential psi (V)");
                        if show_legends {
                            plot_psi =
                                plot_psi.legend(Legend::default().position(Corner::LeftBottom));
                        }
                        plot_psi.show(ui, |plot_ui| {
                            plot_ui.line(
                                Line::new(PlotPoints::new(pts_psi))
                                    .name("Electrostatic Potential psi (V)")
                                    .color(Color32::from_rgb(175, 115, 255))
                                    .width(2.0_f32),
                            );
                            plot_ui.vline(
                                VLine::new(xj_um)
                                    .name("Metallurgical Junction xj")
                                    .color(Color32::from_gray(200))
                                    .style(LineStyle::Dashed { length: 6.0 }),
                            );
                        });
                    });
                });
            }
        }
    }

    /// Renders Space Charge Density rho(x) [C/cm^3]
    fn draw_space_charge_density(
        ui: &mut egui::Ui,
        res: &SimulationResult,
        xj_um: f64,
        show_legends: bool,
    ) {
        let n = res.x_um.len();
        let pts_rho: Vec<[f64; 2]> = (0..n).map(|i| [res.x_um[i], res.rho_c_cm3[i]]).collect();

        let mut plot = Plot::new("plot_space_charge")
            .x_axis_label("Position x (um)")
            .y_axis_label("Space Charge Density rho (C/cm^3)");

        if show_legends {
            plot = plot.legend(Legend::default().position(Corner::LeftBottom));
        }

        plot.show(ui, |plot_ui| {
            plot_ui.line(
                Line::new(PlotPoints::new(pts_rho))
                    .name("Charge Density rho (C/cm^3)")
                    .color(Color32::from_rgb(255, 80, 180))
                    .width(2.0_f32),
            );
            plot_ui.hline(
                HLine::new(0.0)
                    .color(Color32::from_gray(120))
                    .style(LineStyle::Dashed { length: 5.0 }),
            );
            plot_ui.vline(
                VLine::new(xj_um)
                    .name("Metallurgical Junction xj")
                    .color(Color32::from_gray(200))
                    .style(LineStyle::Dashed { length: 6.0 }),
            );
        });
    }

    /// Renders I-V Characteristic Curve: log10(|J|) vs V and Ideality Factor eta(V)
    fn draw_iv_characteristic(
        ui: &mut egui::Ui,
        res: &SimulationResult,
        v_bias: f64,
        show_legends: bool,
    ) {
        let iv = &res.iv_curve;
        let n = iv.v_v.len();

        let pts_jtot: Vec<[f64; 2]> = (0..n)
            .map(|i| [iv.v_v[i], iv.j_total_a_cm2[i].log10()])
            .collect();
        let pts_jdiff: Vec<[f64; 2]> = (0..n)
            .map(|i| [iv.v_v[i], iv.j_diff_a_cm2[i].log10()])
            .collect();
        let pts_jsrh: Vec<[f64; 2]> = (0..n)
            .map(|i| [iv.v_v[i], iv.j_srh_a_cm2[i].log10()])
            .collect();

        let pts_eta: Vec<[f64; 2]> = (0..n)
            .filter(|&i| iv.v_v[i] >= 0.05)
            .map(|i| [iv.v_v[i], iv.ideality[i]])
            .collect();

        let avail_h = ui.available_height();
        let top_h = (avail_h - 10.0) * 0.65;
        let bot_h = (avail_h - 10.0) * 0.35;

        ui.vertical(|ui| {
            // Upper plot: Semi-log I-V Characteristic
            ui.allocate_ui(Vec2::new(ui.available_width(), top_h), |ui| {
                let mut plot = Plot::new("plot_iv_semilog")
                    .x_axis_label("Applied Bias V (V)")
                    .y_axis_label("log10(|J| [A/cm^2])");

                if show_legends {
                    plot = plot.legend(Legend::default().position(Corner::LeftBottom));
                }

                plot.show(ui, |plot_ui| {
                    plot_ui.line(
                        Line::new(PlotPoints::new(pts_jtot))
                            .name("Total Current |J_total|")
                            .color(Color32::from_rgb(255, 235, 100))
                            .width(2.5_f32),
                    );
                    plot_ui.line(
                        Line::new(PlotPoints::new(pts_jdiff))
                            .name("Shockley Diffusion |J_diff|")
                            .color(Color32::from_rgb(0, 200, 255))
                            .width(1.8_f32),
                    );
                    plot_ui.line(
                        Line::new(PlotPoints::new(pts_jsrh))
                            .name("SRH Recombination |J_SRH|")
                            .color(Color32::from_rgb(255, 130, 70))
                            .style(LineStyle::Dashed { length: 6.0 })
                            .width(1.8_f32),
                    );
                    plot_ui.vline(
                        VLine::new(v_bias)
                            .name("Operating Point V_bias")
                            .color(Color32::from_rgb(100, 255, 120))
                            .style(LineStyle::Dashed { length: 5.0 }),
                    );
                });
            });

            ui.add_space(8.0);

            // Lower plot: Ideality Factor eta(V)
            ui.allocate_ui(Vec2::new(ui.available_width(), bot_h), |ui| {
                let mut plot_eta = Plot::new("plot_ideality_factor")
                    .x_axis_label("Forward Bias V (V)")
                    .y_axis_label("Ideality Factor eta");

                if show_legends {
                    plot_eta = plot_eta.legend(Legend::default().position(Corner::LeftBottom));
                }

                plot_eta.show(ui, |plot_ui| {
                    plot_ui.line(
                        Line::new(PlotPoints::new(pts_eta))
                            .name("Ideality Factor eta(V)")
                            .color(Color32::from_rgb(180, 130, 255))
                            .width(2.0_f32),
                    );
                    plot_ui.hline(
                        HLine::new(1.0)
                            .name("Ideal Diffusion Limit (eta = 1)")
                            .color(Color32::from_rgb(0, 200, 255))
                            .style(LineStyle::Dashed { length: 5.0 }),
                    );
                    plot_ui.hline(
                        HLine::new(2.0)
                            .name("SRH Limit (eta = 2)")
                            .color(Color32::from_rgb(255, 130, 70))
                            .style(LineStyle::Dashed { length: 5.0 }),
                    );
                    plot_ui.vline(
                        VLine::new(v_bias)
                            .name("Operating Point V_bias")
                            .color(Color32::from_rgb(100, 255, 120))
                            .style(LineStyle::Dashed { length: 5.0 }),
                    );
                });
            });
        });
    }
}
