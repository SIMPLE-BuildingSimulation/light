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
// use geometry3d::{Point3D};
use geometry3d::intersection::SurfaceSide;
use rendering::ray::Ray;
use rendering::interaction::Interaction;
use solar::ReinhartSky;
use matrix::Matrix;
use rendering::rand::*;
use rendering::samplers::HorizontalCosineWeightedHemisphereSampler;
use geometry3d::Vector3D;


#[cfg(feature = "parallel")]
use rayon::prelude::*;


pub fn calc_dc(rays: &[Ray3D], scene: &Scene, mf: usize)-> Matrix {
    
    // let counter = std::sync::Arc::new(std::sync::Mutex::new(0));

    // Initialize DC Factory
    let mut factory = DCFactory::new(mf);
    factory.max_depth = 3;
    factory.n_ambient_samples = 3000;
    
    
    // Initialize matrix
    let n_bins = factory.reinhart.n_bins;    

    // Process... This can be in parallel, or not.
    #[cfg(not(feature = "parallel"))]
    let aux_iter = rays.iter();
    #[cfg(feature = "parallel")]
    let aux_iter = rays.par_iter();
    // Iterate the rays
    let dcs : Vec<Matrix> = aux_iter.map(|ray|-> Matrix {
        
        // let normal = ray.direction;        
        let origin = ray.origin;        
        
        // Run each spawned ray in parallel or series, depending on 
        // the compilation options
        
        let aux_iter = HorizontalCosineWeightedHemisphereSampler::new(factory.n_ambient_samples);
        #[cfg(feature = "parallel")]
        let aux_iter = {
            let aux : Vec<Vector3D>= aux_iter.map(|v|v).collect();
            aux.into_par_iter()
        };
        
        
        
        // Iterate primary rays           
        let ray_contributions : Vec<Matrix> = aux_iter.map(| new_ray_dir: Vector3D| -> Matrix {
            let mut this_ret = Matrix::new(0.0, 1, n_bins);

            debug_assert!((1.-new_ray_dir.length()).abs() < 0.0000001);
            let new_ray = Ray{
                time: 0.,
                geometry: Ray3D {
                    direction : new_ray_dir,
                    origin,
                }
            };

            
            let mut rng = get_rng();
            // let current_weight = cos_theta;
            factory.trace_ray(scene, &new_ray, 0, PI, factory.n_ambient_samples, &mut this_ret, &mut rng);                                                
            
            
            // let mut c = counter.lock().unwrap();
            // *c += 1;
            // let nrays = rays.len() * factory.n_ambient_samples;
            // let perc = (100. *  *c as Float/ nrays  as Float).round() as usize;            
            // eprintln!("Ray {} of {} ({}%) done...", c, nrays, perc);
            
            this_ret
        }).collect();// End of iterating primary rays

        let mut ret = Matrix::new(0.0, 1, n_bins);     
        ray_contributions.iter().for_each(|v| {
            ret.add_to_this(&v).unwrap();
        });
        ret
    }).collect(); // End of iterating rays
    

    // Write down the results
    let mut ret = Matrix::new(0.0, rays.len(), n_bins);    
    for (sensor_index, contribution) in dcs.iter().enumerate(){
        // add contribution                 
        for patch_index in 0..n_bins{
            let v = contribution.get(0, patch_index).unwrap();
            ret.set(sensor_index, patch_index, v ).unwrap();
        }                                             
    }

    ret
}

/// A structure meant to calculate DC matrices
/// for Climate Daylight Simulations.
pub struct DCFactory {
    pub reinhart: ReinhartSky,
    pub max_depth: usize,    
    pub n_ambient_samples: usize,
    pub limit_weight: Float,
    // pub limit_reflections: usize,
}

impl Default for DCFactory{
    fn default()->Self{
        Self{
            reinhart: ReinhartSky::new(1),
            max_depth: 0,            
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
     fn trace_ray(&self, scene: &Scene, ray: &Ray, current_depth: usize, current_value: Float,  denom_samples: usize, contribution: &mut Matrix, rng: &mut RandGen){
        // Limit bounces        
        if current_depth > self.max_depth {            
            return 
        }
        
        let one_over_samples = 1./ self.n_ambient_samples as Float;        
        // If hits an object
        if let Some((t, interaction)) = scene.cast_ray(ray) {            
            if let Interaction::Surface(data) = &interaction{
                let object = &scene.objects[data.prim_index];
                // get the normal... can be textured.           
                                
                
                let normal = data.geometry_shading.normal.get_normalized();
                let e1 = data.geometry_shading.dpdu.get_normalized();
                let e2 = e1.cross(normal).get_normalized();
                debug_assert!((1.0 - normal.length()).abs() < 0.000001);
                
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
                
                let intersection_pt = ray.geometry.project(t);
                let ray_dir = ray.geometry.direction;
            
                // for now, emmiting materials don't reflect
                if !material.emits_direct_light() {
                    
                    // Run each spawned ray                    
                    
                    /* Adapted From Radiance's samp_hemi() at src/rt/ambcomp.c */
                    let mut wt = current_value;
                    
                    let d = 0.8* current_value * current_value * one_over_samples / self.limit_weight;
                    if wt > d {
                        wt = d;
                    }
                    let mut n = ((self.n_ambient_samples as Float * wt).sqrt() + 0.5).round() as usize;                    
                    if n < 1 {
                        n = 1;
                    }
                    /* End of Adapted Radiance's code*/
                    
                    (0..n).for_each(|_| {
                    
                        let (new_ray_dir, _material_pdf, _is_specular) = material.sample_bsdf(normal, e1, e2, ray_dir, rng);                            
                        

                        debug_assert!((1. as Float-new_ray_dir.length()).abs() < 1e-5, "Length is {}", new_ray_dir.length());
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
                        
                        let new_value = current_value * cos_theta * refl */*material_pdf * one_over_samples * */ 1.5;

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
                        }
                        self.trace_ray(scene, &new_ray, current_depth + 1, new_value, denom_samples * n, contribution, rng);                            
                        
                    });// End the foreach spawned ray
                }
            } else {
                unreachable!();
            }
                        
        } else {        

            let bin_n = self.reinhart.dir_to_bin(ray.geometry.direction);
            let li = 1.;
            let old_value = contribution.get(0, bin_n).unwrap();
            contribution.set(0,bin_n, old_value + li * current_value / denom_samples as Float).unwrap();
            
            
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
