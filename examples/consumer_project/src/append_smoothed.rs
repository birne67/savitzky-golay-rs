/* This program compares signal data and an applied Savitzky-Golay filter 
    in python and rust.
   It reads a CSV file with at least two columns: "Raw" and "Python".
   It applies the Savitzky-Golay filter to the "Raw" data and appends the smoothed data
   as a new column "smoothed" to the same CSV file.
   Finally, it visualizes the original, python, and smoothed data using egui/eframe.
*/

use savitzky_golay::{read_csv_column, BoundaryMode, SavitzkyGolayFilter};
use std::error::Error;

// Add eframe/egui imports
use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};
use egui::Stroke;
use egui::Color32;

fn main() -> Result<(), Box<dyn Error>> {

    // Define file paths and column names
    let input_path = "data.csv";
    let output_path = "data.csv"; // Overwrite the same file
    let column_name = "Raw"; // Change as needed
    let python_column_name = "Python"; // Change as needed
    let smoothed_col_name = "smoothed";

    // Read the whole CSV into memory
    let mut rdr = csv::Reader::from_path(input_path)?;
    println!("Debug: Starting to read CSV file.");
    let headers = rdr.headers()?.clone();
    println!("Debug: Read headers: {:?}", headers);
    let mut records: Vec<csv::StringRecord> = rdr.records().map(|r| r.unwrap()).collect();
    println!("Debug: Read {} records from {}", records.len(), input_path);

    // Debug: Print all available column names
    for header in headers.iter() {
        println!("Available column: {}", header);
    }

    // Read the raw signal data column as f64
    let data = read_csv_column(input_path, column_name)?;
    println!(
        "Read {} data points from column '{}'",
        data.len(),
        column_name
    );

    //Read the python column as f64 for comparison
    let python_data = read_csv_column(input_path, python_column_name)?;
    println!(
        "Read {} data points from column '{}'",
        python_data.len(),
        python_column_name
    );


    // Apply Savitzky-Golay filter to the signal data
    // 41 is the window size (must be odd), 5 is the polynomial order
    // nearest boundary handling
    // like used in the python example
    
    let mut nearest_filter =
        SavitzkyGolayFilter::new(41, 5)?.with_boundary_mode(BoundaryMode::Nearest);
    let smoothed = nearest_filter.apply(&data);

    // Add the new column to each record
    if headers.iter().any(|h| h == smoothed_col_name) {
        println!("Warning: Column '{}' already exists. No additional column will be added.", smoothed_col_name);
    } else {
        
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
    }  

    // Visualization with egui/eframe
    // Set dark mode visuals
    let native_options = eframe::NativeOptions::default();
    // Dark mode is set in the App's update method.

    // Launch egui plot
    let original = data.clone();
    let smoothed_vec = smoothed.clone();
    let _ = eframe::run_native(
        "Savitzky-Golay Plot",
        native_options,
        Box::new(|_cc| {
            Ok(Box::new(PlotApp {
                original,
                python: python_data,
                smoothed: smoothed_vec,
            }))
        }),
    );
    Ok(())
}

// Define the eframe application structure
struct PlotApp {
    original: Vec<f64>,
    python: Vec<f64>,
    smoothed: Vec<f64>,
}

// Implement the eframe application
impl eframe::App for PlotApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        ctx.set_visuals(egui::Visuals::dark());
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("Original, Python, and Calculated Data");

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

            // The python data points for verification of the calculation
            let points_python_vec: Vec<[f64; 2]> = self
                .python
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

            // The data points with the here applied Savitzky-Golay filter
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

            let has_data = !points_original_vec.is_empty() || !points_smoothed_vec.is_empty();
            if !has_data {
                ui.label("No valid data to plot. Check your input CSV and filter settings.");
                return;
            }

            // Determine plot bounds
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

            // Collect lengths of data series
            let original_len = points_original_vec.len();
            let python_len = points_python_vec.len();
            let smoothed_len = points_smoothed_vec.len();

            // Create the plot
            Plot::new("plot")
                .legend(egui_plot::Legend::default())
                .x_axis_label("Index")
                .y_axis_label("Value")
                .include_x(min_x)
                .include_x(max_x)
                .include_y(min_y)
                .include_y(max_y)
                .show(ui, |plot_ui| {
                    // Plot the original data points
                    if original_len > 0 {
                        plot_ui.line(Line::new(
                            "Original",
                            PlotPoints::from(points_original_vec.clone()),
                        )
                        .stroke(Stroke::new(2.0, Color32::GRAY))
                        .style(egui_plot::LineStyle::Solid));
                    }
                    // Plot the python data points
                    if python_len > 0 {
                        plot_ui.line(Line::new(
                            "Python",
                            PlotPoints::from(points_python_vec.clone()),
                        )
                        .stroke(Stroke::new(2.0, Color32::RED))
                        .style(egui_plot::LineStyle::Solid));
                    }
                    // Plot the smoothed data points
                    if smoothed_len > 0 {
                        plot_ui.line(Line::new(
                            "Smoothed",
                            PlotPoints::from(points_smoothed_vec.clone()),
                        )
                        .color(Color32::BLUE)
                        .style(egui_plot::LineStyle::Dotted { spacing: 5.0 }));
                    }
                });

            });
    }
}