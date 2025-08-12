use nalgebra::{DMatrix};
use ndarray:: Array3; 
use std::f64::consts::PI;
use plotters::prelude::*;
use ndarray::s;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // From yees paper 
    let alpha: f64 = 0.025;

    // Resolution Factor
    let factor: f64=5.0;

    // Grid Partitions
    let dx: f64 = alpha / ((factor as f64)*8.0);
    let dy: f64 = alpha / ((factor as f64)*8.0);
    let dtau: f64 = alpha / ((factor as f64) *16.0);

    // Speed of light 
    let cee=3e8;

    // E&M Constants
    let Zee: f64=376.7;
    //Time Partition
    let dt = (alpha)/(factor*16.0*cee);

    // Grid Dimensions (same as Yees)
    let nx: usize = 81;
    let ny: usize = 98;
    let nt: usize = 200;
    let xaxis: Vec<f64> = (0..nx).map(|i|i as f64).collect();

    // Hx and Hy are not plotted
    // They will be 2-dimension and not stored 
    let Hxdata: Vec<f64> = vec![0.0; ny * nx]; 
    let mut Hx: DMatrix<f64> = DMatrix::from_row_slice(ny,nx, &Hxdata);

    let Hydata: Vec<f64> = vec![0.0; ny*nx];
    let mut Hy : DMatrix<f64>= DMatrix::from_row_slice(ny,nx, &Hydata);

    // Create 3-dimensional array for Ez 
    // Will plot at varying locations on x and at different time, t
    let mut Ez: Array3<f64> = Array3::<f64>::zeros((nt,ny,nx));

    // Tuning parameters for incidental wave
    // Location
    let gee: f64=1.5;

    // Width
    let hhh: f64 =0.5;

    //Creating the initial conditions at t=0
    for i in 0..nx{
        let x: f64 = (i as f64) * dx;
        let aux: f64 = x - gee * alpha;
        for j in 0..ny{
            if (0.0 <= aux) && (aux <= (hhh * alpha)){
                Ez[[0,j,i]] = ((aux * PI)/(hhh*alpha)).sin();
            }
        }
    }

    // Obstacle boundary conditions
    for i in 17..49{
        for j in 33..65{
                Ez[[0,j,i]] = 0.0;
            
        }
    }

    for n in 0..(nt-1){
        let time =(n as f64)*dt;
        //Applying the incident wave first 
        for i in 0..nx{
            let x = (i as f64) * dx;
            let aux = x- gee * alpha + cee * time;
            for j in 0..ny {
                if 0.0 <= aux && aux <= (hhh * alpha) {
                    Ez[[n, j, i]] = ((aux * PI) / (hhh * alpha)).sin();
                }
            }
        }
        
        // Ensuring Obstacle Boundary conditions are enforced
        for i in 17..49{
            for j in 33..65{
                    Ez[[n,j,i]] = 0.0;     
            }
        } 
        //Applying HX individually due to grid indexing
        for i in 0..nx{
            for j in 0..ny-1{
                Hx[(j,i)] -= (1.0/Zee)*(dtau/dy)*(Ez[[n,(j+1),i]]-Ez[[n,j,i]]);
            }
        }
        
        //Applying HY individually due to grid indexing
        for i in 0..(nx-1){
            for j in 0..ny{
                Hy[(j,i)] +=(1.0/Zee)*(dtau/dx)*(Ez[[n,j,(i+1)]]-Ez[[n,j,i]])
            }
        }

        //Calculating and storing the value for Ez to be plotted 
        for i in 1..(nx-1){
            for j in 1..(ny-1){
                if j+1<ny{
                    Ez[[(n+1),j,i]] = Ez[[n,j,i]]+Zee*(dtau/dx)*(Hy[(j,i)]-Hy[(j,(i-1))])-Zee*(dtau/dy)*(Hx[(j,i)]-Hx[((j-1),i)]);
                }
            }
        }
    
        // Edge Boundary Conditions 
        for j in 0..ny {
            Ez[[n+1,j, 0]] = 0.0;
            Ez[[n, j, (nx-1)]] = 0.0;
            }

        for i in 0..nx{
            Ez[[n+1,0, i]] = 0.0;
            Ez[[n, (ny-1), i]] = 0.0;
            }    

        // Obstacle Boundary conditions
        for i in 17..49{
            for j in 33..65{
                    Ez[[n+1,j,i]] = 0.0;                
            }
        } 
        
    }
    
    // Plotting
    let n: i32 = 195;  
    let j: i32 = 50; // fixed j value
    let row = Ez.slice(s![n, j, ..]).to_owned(); // shape: (nx,)

    let (min_date, max_date) = xaxis
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), d| {
            (min.min(*d), max.max(*d))
        });

    let (min_alt, max_alt) = row
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), a| {
            (min.min(*a), max.max(*a))
        });

        // Set up drawing backend
    let root = BitMapBackend::new("No_Obstacle.png", (700,500)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Ez versus Position", ("sans-serif", 20))
        .margin(30)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(min_date..max_date, min_alt..max_alt)?;

    chart.configure_mesh()
        .x_desc("Position (x)")
        .y_desc("Ez")
        .draw()?;

    let cleaned_data: Vec<(f64, f64)> = xaxis
        .iter()
        .zip(row.iter())
        .map(|(&x, &y)| (x, y))
        .collect();

    chart.draw_series(LineSeries::new(cleaned_data, &BLUE))?;

    Ok(())   
}