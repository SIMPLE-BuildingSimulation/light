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
use rand::prelude::*;

pub fn calc_dc(rays: &Vec<Ray3D>, scene: &Scene, mf: usize)-> Matrix {
    // Initialize DC Factory
    let mut factory = DCFactory::new(mf);
    factory.max_depth = 1;
    factory.n_ambient_samples = 1000;
    factory.limit_weight=0.0;

    // Initialize matrix
    let n_bins = factory.reinhart.n_bins;
    let mut ret = Matrix::new(0.0, rays.len(), n_bins);

    // Process.
    for (sensor_index, ray) in rays.iter().enumerate(){
        let mut rng = rand::thread_rng();
        let normal = ray.direction;
        let e2 = normal.get_perpendicular().unwrap();
        let e1 = e2.cross(normal);
        let origin = ray.origin;

        let mut spectrum : Vec<Float> = vec![0.0; n_bins];

        for _ in 0..factory.n_ambient_samples {
            // Choose a direction.            
            let new_ray_dir = rendering::samplers::cosine_weighted_sample_hemisphere(&mut rng, e1, e2, normal);
            debug_assert!((1.-new_ray_dir.length()).abs() < 0.0000001);
            let new_ray = Ray{
                time: 0.,
                geometry: Ray3D {
                    direction : new_ray_dir,
                    origin,
                }
            };

            let cos_theta = (normal * new_ray_dir).abs();
            let current_prob = cos_theta;
            let current_weight = cos_theta;
            factory.trace_ray(&mut rng, scene, &new_ray, 0, current_weight, current_prob, &mut spectrum);                        

            
            // // add contribution                                      
            let one_over_samples = 1./ factory.n_ambient_samples as Float;
            for (patch_index, v) in spectrum.iter().enumerate(){
                ret.set(sensor_index, patch_index, *v / one_over_samples).unwrap();
            }
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
    pub limit_reflections: usize,
}

impl Default for DCFactory{
    fn default()->Self{
        Self{
            reinhart: ReinhartSky::new(1),
            max_depth: 2,
            n_shadow_samples: 900,
            n_ambient_samples: 10,

            limit_weight: 1e-5,
            limit_reflections: 0,
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
     fn trace_ray(&self, rng: &mut ThreadRng, scene: &Scene, ray: &Ray, current_depth: usize, current_value: Float, current_prob: Float, spectrum: &mut Vec<Float>){
        
        let n_samples = self.n_ambient_samples as Float;
        let one_over_samples_sq = 1./n_samples/n_samples;

        // Limit bounces        
        if current_depth > self.max_depth {
            return 
        }

        
        
        
        // If hits an object
        if let Some((t, interaction)) = scene.cast_ray(ray) {

            let object = interaction.object();
            if let Interaction::Surface(data) = &interaction{
                // get the normal... can be textured.           
                let normal = data.normal();

                debug_assert!((1.0 - normal.length()).abs() < 0.000001);

                // let object = interaction.object();

                let material = match data.geometry_shading.side {
                    SurfaceSide::Front => {
                        &scene.materials[object.front_material_index]
                    },
                    SurfaceSide::Back =>{
                        &scene.materials[object.back_material_index]
                    }                        
                };
                
                let intersection_pt = ray.geometry.project(t);
                
                let ray_dir = ray.geometry.direction;
            
                
                // for now, emmiting materials don't reflect
                if !material.emits_direct_light() {
                    
                    
                    for _ in 0..self.n_ambient_samples {
                        // Choose a direction.
                        let new_ray_dir = material.sample_bsdf(rng, ray_dir, data.geometry_shading);
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
                        
                        let material_pdf = material.bsdf(ray_dir, normal, new_ray_dir);
                        let new_value = current_value * cos_theta * refl * material_pdf;

                        // Check reflection limits... as described in RTRACE's man
                        if new_value < self.limit_weight && self.limit_reflections > 0 {
                            return;
                        } else {
                            // russian roulette
                            let q : Float = rng.gen();
                            if q > new_value/self.limit_weight {
                                return;
                            }
                        }

                        let new_prob = current_prob * material_pdf;

                        self.trace_ray(rng, scene, &new_ray, current_depth + 1, new_value, new_prob, spectrum);                            
                    }                        
                }
            }else{
                unreachable!();
            }
                        
        } else {        

            
            let bin_n = self.reinhart.dir_to_bin(ray.geometry.direction);
            
            spectrum[bin_n] += PI * one_over_samples_sq * current_value /current_prob ;
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
