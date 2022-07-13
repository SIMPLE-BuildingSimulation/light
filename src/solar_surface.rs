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
use std::path::Path;
use std::rc::Rc;

use crate::Float;
use matrix::Matrix;
use rendering::{colour_matrix::*, DCFactory, Scene};

use simple_model::{Fenestration, SimulationStateElement, SimulationStateHeader, Surface};

use geometry3d::{
    Point3D,
    Polygon3D,
    Ray3D,
    Triangulation3D,
    Vector3D,
    // Loop3D,
};
use rendering::primitive_samplers::sample_triangle_surface;
use rendering::rand::*;

fn get_sampler(triangles_areas: Vec<Float>) -> impl Fn(&mut RandGen) -> usize {
    let total_area: Float = triangles_areas.iter().sum();

    move |rng: &mut RandGen| -> usize {
        let mut r: Float = rng.gen();
        r *= total_area;
        let mut acc = 0.0;
        for (i, area) in triangles_areas.iter().enumerate() {
            acc += area;
            if r <= acc {
                return i;
            }
        }
        unreachable!();
    }
}

/// Structure that can help calculate solar radiation
///
/// It contains the normal of the original Surface and the points
/// randomly sampled in each surface.
pub struct SolarSurface {
    points: Vec<Point3D>,
    pub normal: Vector3D,
    // nrays: usize,
}

impl SolarSurface {
    /// Offset for the starting point of the rays.
    const DELTA: Float = 0.001;

    pub fn new(nrays: usize, polygon: &Polygon3D) -> Self {
        // Get polygon
        let normal = polygon.normal();

        // Triangulate the polygon
        let triangles = Triangulation3D::from_polygon(polygon)
            .unwrap()
            .get_trilist();
        let triangles_areas: Vec<Float> = triangles.iter().map(|t| t.area()).collect();

        // Build a triangle sampler
        let sampler = get_sampler(triangles_areas);
        let mut rng = get_rng();

        // sample points
        let points: Vec<Point3D> = (0..nrays)
            .map(|_| {
                // choose the triangle
                let i = sampler(&mut rng);
                // choose a point in the triangle
                sample_triangle_surface(&triangles[i], &mut rng)
            })
            .collect();

        // return
        Self {
            normal,
            points,
            // nrays,
        }
    }

    /// Receives a list of `SolarSurface` objects as well as the `Scene` containing them and
    /// calculates the back Daylight Coefficient Matrix that can be used for
    /// estimating the incident solar radiation in W/m2. The options for this calculation are
    /// contained in the `DCFactory` used as input.
    fn calc_back_solar_dc_matrix(
        list: &[SolarSurface],
        scene: &Scene,
        dc_factory: &DCFactory,
    ) -> Matrix {
        if list.is_empty() {
            return Matrix::empty();
        }

        // Then the back
        let back_dcs: Vec<Matrix> = list
            .iter()
            .map(|s| s.back_solar_irradiance(scene, dc_factory))
            .collect();
        let mut back = back_dcs[0].clone();
        for dc in back_dcs.iter().skip(1) {
            back.concat_rows(dc).unwrap();
        }
        // return
        back
    }

    /// Receives a list of `SolarSurface` objects as well as the `Scene` containing them and
    /// calculates the front Daylight Coefficient Matrix that can be used for
    /// estimating the incident solar radiation in W/m2. The options for this calculation are
    /// contained in the `DCFactory` used as input.
    fn calc_front_solar_dc_matrix(
        list: &[SolarSurface],
        scene: &Scene,
        dc_factory: &DCFactory,
    ) -> Matrix {
        if list.is_empty() {
            return Matrix::empty();
        }

        // Calculate the front
        let front_dcs: Vec<Matrix> = list
            .iter()
            .map(|s| s.front_solar_irradiance(scene, dc_factory))
            .collect();
        let mut front = front_dcs[0].clone();
        for dc in front_dcs.iter().skip(1) {
            front.concat_rows(dc).unwrap();
        }

        // return
        front
    }

    /// Gets the front Daylight Coefficient Matrix that can be used for
    /// estimating the incident solar radiation in W/m2. If the `path` exists,
    /// then it will attempt to read the matrix in that file. If it does not exist,
    /// it will calculate it and write it down there. If no `path` is given, it will
    /// just calculate and not save it anywhere.
    pub fn get_front_solar_dc_matrix(
        list: &[SolarSurface],
        path: Option<&String>,
        scene: &Scene,
        dc_factory: &DCFactory,
    ) -> Result<Matrix, String> {
        let matrix = if let Some(path) = path {
            let path = Path::new(path);
            if path.exists() && path.is_file() {
                // Attempt to read... return error if error
                read_matrix(path)?
            } else {
                // Calculate and write it down
                let m = Self::calc_front_solar_dc_matrix(list, scene, dc_factory);
                save_matrix(&m, path)?;
                m
            }
        } else {
            // Just calculate
            Self::calc_front_solar_dc_matrix(list, scene, dc_factory)
        };
        // return
        Ok(matrix)
    }

    /// Gets the back Daylight Coefficient Matrix that can be used for
    /// estimating the incident solar radiation in W/m2. If the `path` exists,
    /// then it will attempt to read the matrix in that file. If it does not exist,
    /// it will calculate it and write it down there. If no `path` is given, it will
    /// just calculate and not save it anywhere.
    pub fn get_back_solar_dc_matrix(
        list: &[SolarSurface],
        path: Option<&String>,
        scene: &Scene,
        dc_factory: &DCFactory,
    ) -> Result<Matrix, String> {
        let matrix = if let Some(path) = path {
            let path = Path::new(path);
            if path.exists() && path.is_file() {
                // Attempt to read... return error if error
                read_matrix(path)?
            } else {
                // Calculate and write it down
                let m = Self::calc_back_solar_dc_matrix(list, scene, dc_factory);
                save_matrix(&m, path)?;
                m
            }
        } else {
            // Just calculate
            Self::calc_back_solar_dc_matrix(list, scene, dc_factory)
        };
        // return
        Ok(matrix)
    }

    /// Builds a set of SolarSurfaces from Fenestrations
    ///
    /// Adds the necessary elements to the `SimulationStateHeader`
    pub fn make_fenestrations(
        list: &[Rc<Fenestration>],
        state: &mut SimulationStateHeader,
        n_rays: usize,
    ) -> Vec<SolarSurface> {
        list.iter()
            .enumerate()
            .map(|(i, s)| {
                let i = state.push(
                    SimulationStateElement::FenestrationFrontSolarIrradiance(i),
                    0.0,
                );
                s.set_front_incident_solar_irradiance_index(i);

                let i = state.push(
                    SimulationStateElement::FenestrationBackSolarIrradiance(i),
                    0.0,
                );
                s.set_back_incident_solar_irradiance_index(i);

                let i = state.push(
                    SimulationStateElement::FenestrationFrontIRIrradiance(i),
                    0.0,
                );
                s.set_front_ir_irradiance_index(i);

                let i = state.push(SimulationStateElement::FenestrationBackIRIrradiance(i), 0.0);
                s.set_back_ir_irradiance_index(i);
                // Create
                SolarSurface::new(n_rays, &s.vertices)
            })
            .collect()
    }

    /// Builds a set of SolarSurfaces from Surfaces
    ///
    /// Adds the necessary elements to the `SimulationStateHeader
    pub fn make_surfaces(
        list: &[Rc<Surface>],
        state: &mut SimulationStateHeader,
        n_rays: usize,
    ) -> Vec<SolarSurface> {
        list.iter()
            .enumerate()
            .map(|(i, s)| {
                let i = state.push(SimulationStateElement::SurfaceFrontSolarIrradiance(i), 0.0);
                s.set_front_incident_solar_irradiance_index(i);

                let i = state.push(SimulationStateElement::SurfaceBackSolarIrradiance(i), 0.0);
                s.set_back_incident_solar_irradiance_index(i);

                let i = state.push(SimulationStateElement::SurfaceFrontIRIrradiance(i), 0.0);
                s.set_front_ir_irradiance_index(i);

                let i = state.push(SimulationStateElement::SurfaceBackIRIrradiance(i), 0.0);
                s.set_back_ir_irradiance_index(i);

                // create
                SolarSurface::new(n_rays, &s.vertices)
            })
            .collect()
    }

    /// Gets the front rays of a surface
    pub fn front_rays(&self) -> Vec<Ray3D> {
        self.points
            .iter()
            .map(|p| Ray3D {
                direction: self.normal,
                origin: *p + self.normal * Self::DELTA,
            })
            .collect()
    }

    /// Gets the back rays of a surface
    pub fn back_rays(&self) -> Vec<Ray3D> {
        self.points
            .iter()
            .map(|p| Ray3D {
                direction: self.normal * -1.,
                origin: *p - self.normal * Self::DELTA,
            })
            .collect()
    }

    /// Calculates the Daylight Coefficient matrix for the front of a `SolarSurface`
    pub fn front_solar_irradiance(&self, scene: &Scene, factory: &DCFactory) -> Matrix {
        let front_rays = self.front_rays();
        let dc = factory.calc_dc(&front_rays, scene);
        let dc = colour_matrix_to_radiance(&dc);
        average_matrix(&dc)
    }

    /// Calculates the Daylight Coefficient matrix for the back of a `SolarSurface`
    pub fn back_solar_irradiance(&self, scene: &Scene, factory: &DCFactory) -> Matrix {
        let back_rays = self.back_rays();
        let dc = factory.calc_dc(&back_rays, scene);
        let dc = colour_matrix_to_radiance(&dc);
        average_matrix(&dc)
    }
}
