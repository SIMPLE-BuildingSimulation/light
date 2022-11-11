/*
MIT License
Copyright (c) 2021 Germán Molina
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
use solar::ReinhartSky;
// use rendering::from_radiance::from
use clap::Parser;
use geometry3d::{Point3D, Ray3D, Vector3D};
use rendering::{DCFactory, Scene};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Inputs {
    #[arg(short, long)]
    input: String,
    // #[arg(short, long)]
    // weather: String,
}

fn main() {
    let args = Inputs::parse();

    let input_file = args.input;
    // let _weather_file = "asd"; //matches.value_of("weather").unwrap();

    let mut scene = if input_file.ends_with(".rad") {
        Scene::from_radiance(input_file)
    } else if input_file.ends_with(".simple") || input_file.ends_with(".spl") {
        panic!("Reading SIMPLE models is still not suppoerted")
    } else {
        eprintln!(
            "Don't know how to read file '{}'... only .rad is supported for now",
            input_file
        );
        std::process::exit(1);
    };

    scene.build_accelerator();

    // Setup sensors
    let up = Vector3D::new(0., 0., 1.);
    let rays = vec![
        Ray3D {
            origin: Point3D::new(2., 0.5, 0.8),
            direction: up,
        },
        Ray3D {
            origin: Point3D::new(2., 2.5, 0.8),
            direction: up,
        },
        Ray3D {
            origin: Point3D::new(2., 5.5, 0.8),
            direction: up,
        },
    ];

    eprintln!("Ready to calc!... # Surface = {}", scene.triangles.len());

    let mf = 1;
    let factory = DCFactory {
        max_depth: 1,
        n_ambient_samples: 10000,
        reinhart: ReinhartSky::new(mf),
        ..DCFactory::default()
    };
    let dc_matrix = factory.calc_dc(&rays, &scene);
    println!("Matrix = {}", dc_matrix);
}
