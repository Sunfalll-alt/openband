//! openband: Native Linux GUI Desktop Application for Interactive
//! 1D Semiconductor Device Simulation (p-n junctions and heterojunctions).
//!
//! Uses eframe's glow (OpenGL) backend with native Wayland and X11 display support.

use eframe::egui::ViewportBuilder;
use openband::OpenBandApp;

fn main() -> eframe::Result<()> {
    // Standard Linux environment configuration
    let native_options = eframe::NativeOptions {
        viewport: ViewportBuilder::default()
            .with_inner_size([1440.0, 900.0])
            .with_min_inner_size([1024.0, 680.0])
            .with_title("openband — 1D Semiconductor Device Simulator"),
        ..Default::default()
    };

    eframe::run_native(
        "openband",
        native_options,
        Box::new(|cc| Ok(Box::new(OpenBandApp::new(cc)))),
    )
}
