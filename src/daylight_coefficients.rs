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


use crate::PI;
use crate::Float;
use rendering::scene::{Scene};
use geometry3d::Ray3D;
use geometry3d::intersect_trait::SurfaceSide;
use rendering::ray::Ray;
use rendering::interaction::Interaction;
use solar::ReinhartSky;
use matrix::Matrix;
use std::sync::{Mutex, Arc};
use rendering::rand::*;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

// const ONE_OVER_PI : Float = 1./PI;

pub fn calc_dc(rays: &[Ray3D], scene: &Scene, mf: usize)-> Matrix {
    // Initialize DC Factory
    let mut factory = DCFactory::new(mf);
    factory.max_depth = 3;
    factory.n_ambient_samples = 600;
    factory.limit_weight= 1e-32; // down to about 1e-40?

    let one_over_samples = 1./ factory.n_ambient_samples as Float;
    
    // Initialize matrix
    let n_bins = factory.reinhart.n_bins;
    let mut ret = Matrix::new(0.0, rays.len(), n_bins);

    // Process... This can be in parallel, or not.
    #[cfg(not(feature = "parallel"))]
    let aux_iter = rays.iter();
    #[cfg(feature = "parallel")]
    let aux_iter = rays.par_iter();

    let dcs : Vec<Arc<Mutex<Vec<Float>>>> = aux_iter.map(|ray|-> Arc<Mutex<Vec<Float>>> {
        
        // let normal = ray.direction;        
        let origin = ray.origin;
        
        let spectrum  = Arc::new(Mutex::new(vec![0.0; n_bins]));
        
        // Run each spawned ray in parallel
        #[cfg(not(feature = "parallel"))]
        let aux_iter = 0..factory.n_ambient_samples;
        #[cfg(feature = "parallel")]
        let aux_iter = (0..factory.n_ambient_samples).into_par_iter();
        // let rng_src = get_rng();

        let counter = Arc::new(Mutex::new(0));
        let _ = aux_iter.map(|_|{
                        
            let mut rng = get_rng();//clone_rng(&rng_src);
            // Choose a direction.            
            let new_ray_dir = rendering::samplers::sample_cosine_weighted_horizontal_hemisphere(&mut rng);
            // let cos_theta = (normal * new_ray_dir).abs();
            // let _current_prob = cos_theta*ONE_OVER_PI;

            debug_assert!((1.-new_ray_dir.length()).abs() < 0.0000001);
            let new_ray = Ray{
                time: 0.,
                geometry: Ray3D {
                    direction : new_ray_dir,
                    origin,
                }
            };

            
            // let current_weight = cos_theta;
            factory.trace_ray(scene, &new_ray, 0, PI*one_over_samples/*current_weight*/,  Arc::clone(&spectrum));                                                
            
            // Divide by the probability of the first ray
            // let mut s = spectrum.lock().unwrap();
            // for v in s.iter_mut(){
            //     *v *= one_over_samples;//current_prob;
            // }
            let mut c = counter.lock().unwrap();
            *c += 1;
            
            let perc = (100. *  *c as Float/ factory.n_ambient_samples as Float).round() as usize;            
            eprintln!("Ray {} of {} ({}%) done...", c, factory.n_ambient_samples, perc);
            
            
        }).collect::<()>();
        spectrum
}).collect();
    

    // Write down the results
    
    for (sensor_index, spectrum) in dcs.iter().enumerate(){
        // add contribution            
        let s = &*spectrum.lock().unwrap();                                  
        for (patch_index, v) in s.iter().enumerate(){
            ret.set(sensor_index, patch_index, *v  ).unwrap();
        }
    }

    ret
}

/// A structure meant to calculate DC matrices
/// for Climate Daylight Simulations.
pub struct DCFactory {
    pub reinhart: ReinhartSky,
    pub max_depth: usize,
    pub n_shadow_samples: usize,
    pub n_ambient_samples: usize,

    pub limit_weight: Float,
    // pub limit_reflections: usize,
}

impl Default for DCFactory{
    fn default()->Self{
        Self{
            reinhart: ReinhartSky::new(1),
            max_depth: 0,
            n_shadow_samples: 900,
            n_ambient_samples: 10,

            limit_weight: 1e-5,
            // limit_reflections: 0,
        }
    }
}




impl DCFactory {



    /// Creates a new `DCFactory` with a Reinhart subdivision `mf`
    pub fn new(mf: usize)->Self{
        Self{
            reinhart: ReinhartSky::new(mf),
            .. DCFactory::default()
        }
    }

     /// Recursively traces a ray until it excedes the `max_depth` of the 
     /// `DCFactory` or the ray does not hit anything (i.e., it reaches either
     /// the sky or the ground)
     fn trace_ray(&self, scene: &Scene, ray: &Ray, current_depth: usize, current_value: Float,  spectrum: Arc<Mutex<Vec<Float>>>){
        // eprintln!("A;");
        // Limit bounces        
        if current_depth > self.max_depth {            
            return 
        }
        // eprintln!("B;");
        
        let one_over_samples = 1./ self.n_ambient_samples as Float;
        // eprintln!("C;");
        // If hits an object
        if let Some((t, interaction)) = scene.cast_ray(ray) {
            // eprintln!("D;");
            let object = interaction.object();
            // eprintln!("E;");
            if let Interaction::Surface(data) = &interaction{
                // eprintln!("F;");
                // get the normal... can be textured.           
                let normal = data.normal();
                // eprintln!("G;");
                debug_assert!((1.0 - normal.length()).abs() < 0.000001);
                // let object = interaction.object();
                
                let material = match data.geometry_shading.side {
                    SurfaceSide::Front => {
                        &scene.materials[object.front_material_index]
                    },
                    SurfaceSide::Back =>{
                        &scene.materials[object.back_material_index]
                    },
                    SurfaceSide::NonApplicable => {
                        return;
                    }                   
                };
                // eprintln!("H;");
                
                let intersection_pt = ray.geometry.project(t);
                // eprintln!("I;");
                let ray_dir = ray.geometry.direction;
            
                
                // for now, emmiting materials don't reflect
                if !material.emits_direct_light() {
                    // eprintln!("J;");
                    
                    let bsdf_sampler = material.bsdf_sampler(data.geometry_shading);

                    // Run each spawned ray                    
                    // let rng_src = get_rng();
                    let _ = (0..self.n_ambient_samples).map(|_| {
                        let mut rng = get_rng();//clone_rng(&rng_src);
                        let (new_ray_dir, _material_pdf) = bsdf_sampler(ray_dir, &mut rng);                            

                        debug_assert!((1.-new_ray_dir.length()).abs() < 0.0000001);
                        let new_ray = Ray{
                            time: ray.time,
                            geometry: Ray3D {
                                direction : new_ray_dir,
                                origin: intersection_pt,// + normal * 0.0001, // avoid self shading
                            }
                        };
                        let cos_theta = (normal * new_ray_dir).abs();
                        // WE ARE USING ONLY THE RED COLOR FOR NOW.
                        let refl = material.colour().red;
                        
                        let new_value = current_value * cos_theta * refl */*material_pdf **/ one_over_samples * 1.5;

                        // Check reflection limits... as described in RTRACE's man
                        // if  self.limit_reflections > 0 && new_value < self.limit_weight {
                        //     return;
                        // } else 
                        if self.limit_weight > 0. && new_value < self.limit_weight {
                            
                            // russian roulette
                            let q : Float = rng.gen();
                            if q > new_value/self.limit_weight {                                
                                return;
                            }
                        }else{
                            self.trace_ray(scene, &new_ray, current_depth + 1, new_value,  Arc::clone(&spectrum));                            
                        }
                    }).collect::<()>();                        
                }
            }else{
                unreachable!();
            }
                        
        } else {        

            // eprintln!("BEFORE FINISH 1;");
            let bin_n = self.reinhart.dir_to_bin(ray.geometry.direction);
            // if bin_n == 0 {
            //     eprintln!("Adding {} to bin {}", current_value, bin_n);
            // }
            // eprintln!("BEFORE FINISH 2;");
            let li = 1.;
            let mut s = spectrum.lock().unwrap();
            // eprintln!("BEFORE FINISH 3;");
            // let omega = 2. * PI * one_over_samples;
            // let omega = self.reinhart.bin_solid_angle(bin_n);
            // println!("Oomega = {} | other_omega = {}", omega, self.reinhart.bin_solid_angle(bin_n));
            s[bin_n] +=  li * current_value;// /_current_prob;//*omega ;            
            // eprintln!("BEFORE FINISH 4;");

        }
    }

    
    

    



    
}



#[cfg(test)]
mod tests {
        
    use super::*;
    use geometry3d::{Point3D, Vector3D};
    #[test]
    fn test_calc_dc(){
        // Setup sensors
        let up = Vector3D::new(0., 0., 1.);
        let rays = vec![
            Ray3D{origin: Point3D::new(2., 0.5, 0.8), direction: up },
            Ray3D{origin: Point3D::new(2., 2.5, 0.8), direction: up },
            Ray3D{origin: Point3D::new(2., 5.5, 0.8), direction: up },            
        ];

        // Read scene
        let rad_file = "./test_data/one_surface.rad";
        let scene = Scene::from_radiance(rad_file.to_string());
        eprintln!("Ready to calc!... # Surface = {}", scene.objects.len());
        

        let dc_matrix = calc_dc(&rays, &scene, 1);
        println!("Matrix = {}", dc_matrix);
    }
}