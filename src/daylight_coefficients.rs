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

use crate::colour_matrix::ColourMatrix;
use crate::Float;
use geometry3d::intersection::SurfaceSide;
use geometry3d::Vector3D;
use geometry3d::{Point3D, Ray3D};
use rendering::colour::Spectrum;
use rendering::interaction::Interaction;
use rendering::rand::*;
use rendering::ray::Ray;
use rendering::samplers::HorizontalCosineWeightedHemisphereSampler;
use rendering::scene::Scene;
use solar::ReinhartSky;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// A structure meant to calculate DC matrices
/// for Climate Daylight Simulations.
pub struct DCFactory {
    pub reinhart: ReinhartSky,
    pub max_depth: usize,
    pub n_ambient_samples: usize,
    pub limit_weight: Float,
    // pub limit_reflections: usize,
}

impl Default for DCFactory {
    fn default() -> Self {
        Self {
            reinhart: ReinhartSky::new(1),
            max_depth: 0,
            n_ambient_samples: 10,

            limit_weight: 1e-4,
            // limit_reflections: 0,
        }
    }
}

impl DCFactory {
    pub fn calc_dc(&self, rays: &[Ray3D], scene: &Scene) -> ColourMatrix {
        // Initialize matrix
        let n_bins = self.reinhart.n_bins;

        // Process... This can be in parallel, or not.
        #[cfg(not(feature = "parallel"))]
        let aux_iter = rays.iter();
        #[cfg(feature = "parallel")]
        let aux_iter = rays.par_iter();
        // Iterate the rays
        let dcs: Vec<ColourMatrix> = aux_iter
            .map(|primary_ray| -> ColourMatrix {
                let normal = primary_ray.direction;
                let origin = primary_ray.origin;
                let e2 = normal.get_perpendicular().unwrap();
                let e1 = e2.cross(normal);

                // Run each spawned ray in parallel or series, depending on
                // the compilation options

                let aux_iter =
                    HorizontalCosineWeightedHemisphereSampler::new(self.n_ambient_samples);
                #[cfg(feature = "parallel")]
                let aux_iter = {
                    let aux: Vec<Vector3D> = aux_iter.map(|v| v).collect();
                    aux.into_par_iter()
                };

                // Iterate primary rays
                let ray_contributions: Vec<ColourMatrix> = aux_iter
                    .map(|local_ray_dir: Vector3D| -> ColourMatrix {
                        let (x, y, z) = rendering::samplers::local_to_world(
                            e1,
                            e2,
                            normal,
                            Point3D::new(0., 0., 0.),
                            local_ray_dir.x,
                            local_ray_dir.y,
                            local_ray_dir.z,
                        );
                        let new_ray_dir = Vector3D::new(x, y, z);

                        let mut this_ret = ColourMatrix::new(Spectrum::black(), 1, n_bins);

                        debug_assert!(
                            (1. - new_ray_dir.length()).abs() < 0.0000001,
                            "length is {}",
                            new_ray_dir.length()
                        );

                        let new_ray = Ray {
                            // time: 0.,
                            geometry: Ray3D {
                                direction: new_ray_dir,
                                origin,
                            },
                            refraction_index: 1.,
                        };

                        let mut rng = get_rng();
                        // let current_weight = cos_theta;
                        self.trace_ray(
                            scene,
                            &new_ray,
                            0,
                            Spectrum::gray(1.),
                            self.n_ambient_samples,
                            &mut this_ret,
                            &mut rng,
                        );

                        // let mut c = counter.lock().unwrap();
                        // *c += 1;
                        // let nrays = rays.len() * factory.n_ambient_samples;
                        // let perc = (100. *  *c as Float/ nrays  as Float).round() as usize;
                        // eprintln!("Ray {} of {} ({}%) done...", c, nrays, perc);

                        this_ret
                    })
                    .collect(); // End of iterating primary rays

                let mut ret = ColourMatrix::new(Spectrum::black(), 1, n_bins);
                ray_contributions.iter().for_each(|v| {
                    ret += v;
                });
                ret
            })
            .collect(); // End of iterating rays

        // Write down the results
        let mut ret = ColourMatrix::new(Spectrum::black(), rays.len(), n_bins);
        for (sensor_index, contribution) in dcs.iter().enumerate() {
            // add contribution
            for patch_index in 0..n_bins {
                let v = contribution.get(0, patch_index).unwrap();
                ret.set(sensor_index, patch_index, v).unwrap();
            }
        }

        ret
    }

    /// Recursively traces a ray until it excedes the `max_depth` of the
    /// `DCFactory` or the ray does not hit anything (i.e., it reaches either
    /// the sky or the ground)
    fn trace_ray(
        &self,
        scene: &Scene,
        ray: &Ray,
        current_depth: usize,
        current_value: Spectrum,
        denom_samples: usize,
        contribution: &mut ColourMatrix,
        rng: &mut RandGen,
    ) {
        // Limit bounces
        if current_depth > self.max_depth {
            return;
        }

        // let one_over_samples = 1./ self.n_ambient_samples as Float;
        // If hits an object
        if let Some(interaction) = scene.cast_ray(ray) {
            if let Interaction::Surface(data) = &interaction {
                let object = &scene.objects[data.prim_index];
                // get the normal... can be textured.

                let normal = data.geometry_shading.normal; //.get_normalized();
                let e1 = data.geometry_shading.dpdu.get_normalized();
                let e2 = normal.cross(e1);
                debug_assert!((1.0 - normal.length()).abs() < 0.000001);
                debug_assert!((1.0 - e1.length()).abs() < 0.000001);
                debug_assert!((1.0 - e2.length()).abs() < 0.000001);

                let material = match data.geometry_shading.side {
                    SurfaceSide::Front => &scene.materials[object.front_material_index],
                    SurfaceSide::Back => &scene.materials[object.back_material_index],
                    SurfaceSide::NonApplicable => {
                        return;
                    }
                };

                // let intersection_pt = ray.geometry.project(t);
                let intersection_pt = data.point;
                // let ray_dir = ray.geometry.direction;

                // for now, emmiting materials don't reflect
                if !material.emits_direct_light() {
                    // Run each spawned ray

                    /* Adapted From Radiance's samp_hemi() at src/rt/ambcomp.c */
                    let mut intens = current_value.red;
                    if current_value.green > intens {
                        intens = current_value.green;
                    }
                    if current_value.blue > intens {
                        intens = current_value.blue;
                    }
                    let mut wt = intens;

                    let d = intens * intens * 0.8 /* *one_over_samples*/ / self.limit_weight / denom_samples as Float;
                    if wt > d {
                        wt = d;
                    }
                    let mut n = ((denom_samples as Float * wt).sqrt() + 0.5).round() as usize;
                    if n < 1 {
                        n = 1;
                    }
                    /* End of Adapted Radiance's code*/

                    (0..n).for_each(|_| {
                    
                
                        let (new_ray, _bsdf_value, _is_specular) = material.sample_bsdf(normal, e1, e2, intersection_pt, *ray, rng);                            
                        
                        let new_ray_dir = new_ray.geometry.direction;
                        debug_assert!((1. as Float-new_ray_dir.length()).abs() < 1e-5, "Length is {}", new_ray_dir.length());

                        let cos_theta = (normal * new_ray_dir).abs();
                        let refl = material.colour();

                        // current_value * refl * cos_theta * bsdf_value / pdf == bsdf_value
                        let new_value = current_value * refl * cos_theta  /* * bsdf_value   * one_over_samples * 1.5*/ ;

                        // Check reflection limits... as described in RTRACE's man
                        // if  self.limit_reflections > 0 && new_value < self.limit_weight {
                        //     return;
                        // } 
                        intens = new_value.red;
                        if new_value.green > intens {
                            intens = new_value.green;
                        }
                        if new_value.blue > intens {
                            intens = new_value.blue;
                        }
                        if self.limit_weight > 0. && intens < self.limit_weight {
                            // russian roulette
                            let q : Float = rng.gen();
                            if q > intens/self.limit_weight {
                                return;
                            }
                        }
                        self.trace_ray(scene, &new_ray, current_depth + 1, new_value * 1.5, denom_samples * n, contribution, rng);
                    }); // End the foreach spawned ray
                }
            } else {
                panic!("Found light-emmiting material when calculating Daylight Coefficients");
            }
        } else {
            let bin_n = self.reinhart.dir_to_bin(ray.geometry.direction);

            let li = Spectrum::gray(3.);
            let old_value = contribution.get(0, bin_n).unwrap();

            contribution
                .set(
                    0,
                    bin_n,
                    old_value + li * current_value / denom_samples as Float,
                )
                .unwrap();
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use geometry3d::{Point3D, Vector3D};
    #[test]
    #[ignore]
    fn test_calc_dc() {
        assert!(true);
        // return;
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

        // Read scene
        let rad_file = "./test_data/room.rad";
        let mut scene = Scene::from_radiance(rad_file.to_string());
        scene.build_accelerator();
        eprintln!("Ready to calc!... # Surface = {}", scene.objects.len());

        // Initialize DC Factory
        let factory = DCFactory {
            max_depth: 4,
            n_ambient_samples: 2700,
            reinhart: ReinhartSky::new(1),
            ..DCFactory::default()
        };

        let dc_matrix = factory.calc_dc(&rays, &scene);
        let dc_matrix = crate::colour_matrix::colour_matrix_to_luminance(&dc_matrix);
        crate::colour_matrix::save_matrix(
            &dc_matrix,
            std::path::Path::new("./test_data/full_dc_rust.mtx"),
        )
        .unwrap();
    }
}
