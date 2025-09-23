use savitzky_golay::{read_csv_column, BoundaryMode, SavitzkyGolayFilter};
use std::error::Error;

// Add eframe/egui imports
use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};

fn main() -> Result<(), Box<dyn Error>> {
    let input_path = "data.csv";
    let output_path = "data.csv"; // Overwrite the same file
    let column_name = "Raw"; // Change as needed
    let smoothed_col_name = "smoothed";

    // Read the whole CSV into memory
    let mut rdr = csv::Reader::from_path(input_path)?;
    let headers = rdr.headers()?.clone();
    let mut records: Vec<csv::StringRecord> = rdr.records().map(|r| r.unwrap()).collect();

    println!("Read {} records from {}", records.len(), input_path);

    // Debug: Print all available column names
    println!(
        "Available columns: {:?}",
        headers.iter().collect::<Vec<_>>()
    );

    // Read the target column as f64
    let data = read_csv_column(input_path, column_name)?;
    println!(
        "Read {} data points from column '{}'",
        data.len(),
        column_name
    );

    let mut nearest_filter =
        SavitzkyGolayFilter::new(41, 5)?.with_boundary_mode(BoundaryMode::Nearest);
    let smoothed = nearest_filter.apply(&data);

    // Add the new column to each record
    for (i, record) in records.iter_mut().enumerate() {
        // If smoothed value exists, append it; else, append empty
        let value = smoothed
            .get(i)
            .map(|v: &f64| v.to_string())
            .unwrap_or_default();
        record.push_field(&value);
    }

    // Write out the new CSV
    let mut wtr = csv::Writer::from_path(output_path)?;
    let mut new_headers = headers.clone();
    new_headers.push_field(smoothed_col_name);
    wtr.write_record(&new_headers)?;
    for record in records {
        wtr.write_record(&record)?;
    }
    wtr.flush()?;
    println!("Smoothed column appended to {}!", output_path);

    // Set dark mode visuals
    let mut native_options = eframe::NativeOptions::default();
    native_options.default_theme = eframe::Theme::Dark;

    // Launch egui plot
    let original = data.clone();
    let smoothed_vec = smoothed.clone();
    eframe::run_native(
        "Savitzky-Golay Plot",
        native_options,
        Box::new(|_cc| {
            Box::new(PlotApp {
                original,
                smoothed: smoothed_vec,
            })
        }),
    );
    Ok(())
}

struct PlotApp {
    original: Vec<f64>,
    smoothed: Vec<f64>,
}

impl eframe::App for PlotApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        ctx.set_visuals(egui::Visuals::dark());
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("Original and Smoothed Data");

            // Filter out NaN/inf values for plotting
            let points_original_vec: Vec<[f64; 2]> = self
                .original
                .iter()
                .enumerate()
                .filter_map(|(i, &y)| {
                    if y.is_finite() {
                        Some([i as f64, y])
                    } else {
                        None
                    }
                })
                .collect();
            let points_smoothed_vec: Vec<[f64; 2]> = self
                .smoothed
                .iter()
                .enumerate()
                .filter_map(|(i, &y)| {
                    if y.is_finite() {
                        Some([i as f64, y])
                    } else {
                        None
                    }
                })
                .collect();

            let has_data = points_original_vec.len() > 0 || points_smoothed_vec.len() > 0;
            if !has_data {
                ui.label("No valid data to plot. Check your input CSV and filter settings.");
                return;
            }

            let min_x = 0.0;
            let max_x = self.original.len().max(self.smoothed.len()) as f64;
            let min_y = points_original_vec
                .iter()
                .chain(points_smoothed_vec.iter())
                .map(|p| p[1])
                .fold(f64::INFINITY, f64::min);
            let max_y = points_original_vec
                .iter()
                .chain(points_smoothed_vec.iter())
                .map(|p| p[1])
                .fold(f64::NEG_INFINITY, f64::max);

            let original_len = points_original_vec.len();
            let smoothed_len = points_smoothed_vec.len();

            let points_original: PlotPoints = points_original_vec.into();
            let points_smoothed: PlotPoints = points_smoothed_vec.into();

            Plot::new("plot")
                .legend(egui_plot::Legend::default())
                .x_axis_label("Index")
                .y_axis_label("Value")
                .include_x(min_x)
                .include_x(max_x)
                .include_y(min_y)
                .include_y(max_y)
                .show(ui, |plot_ui| {
                    if original_len > 0 {
                        plot_ui.line(Line::new(points_original).name("Original"));
                    }
                    if smoothed_len > 0 {
                        plot_ui.line(Line::new(points_smoothed).name("Smoothed"));
                    }
                });
        });
    }
}
