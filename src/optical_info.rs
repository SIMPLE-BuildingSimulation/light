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
use crate::solar_surface::SolarSurface;
use crate::Float;
use matrix::Matrix;
use rendering::{DCFactory, Scene, Wavelengths};
use serde::{Deserialize, Serialize};
use simple_model::{SimpleModel, SimulationStateHeader, SolarOptions};
use solar::ReinhartSky;

/// A set of view factors as seen by a `ThermalSurface`.
#[derive(Debug, Default, Copy, Clone, Serialize, Deserialize)]
pub struct IRViewFactorSet {
    /// The fraction of the view that corresponds to the sky
    pub sky: Float,

    /// The fraction of the view that corresponds to the ground
    pub ground: Float,

    /// The fraction of the view that corresponds to other objects and
    /// surfaces (they are assumed to be at air temperature)
    pub air: Float,
}

/// Information about the solar radiation and other optical elements
/// of the whole model.
#[derive(Clone, Serialize, Deserialize)]
pub struct OpticalInfo {
    /// The [`IRViewFactorSet`] for the front side of each surface
    pub front_surfaces_view_factors: Vec<IRViewFactorSet>,

    /// The [`IRViewFactorSet`] for the back side of each surface
    pub back_surfaces_view_factors: Vec<IRViewFactorSet>,

    /// The [`IRViewFactorSet`] for the front side of each fenestration
    pub front_fenestrations_view_factors: Vec<IRViewFactorSet>,

    /// The [`IRViewFactorSet`] for the back side of each fenestration
    pub back_fenestrations_view_factors: Vec<IRViewFactorSet>,

    /// The Daylight Coefficients matrix for the front-side of the  surfaces in the scene
    pub front_surfaces_dc: Matrix,

    /// The Daylight Coefficients matrix for the back-side of the  surfaces in the scene
    pub back_surfaces_dc: Matrix,

    /// The Daylight Coefficients matrix for the front-side of the  fenestrations in the scene
    pub front_fenestrations_dc: Matrix,

    /// The Daylight Coefficients matrix for the back-side of the fenestrations in the scene
    pub back_fenestrations_dc: Matrix,
}

impl OpticalInfo {
    /// Calculates the new OpticalInformation
    pub fn new(
        options: &SolarOptions,
        model: &SimpleModel,
        state: &mut SimulationStateHeader,
    ) -> Result<Self, String> {
        // Collect calculation options
        let mf = *options.solar_sky_discretization().unwrap();
        let n_solar_rays = *options.n_solar_irradiance_points().unwrap();

        // Create Surfaces and Fenestrations
        let surfaces = SolarSurface::make_surfaces(&model.surfaces, state, n_solar_rays)?;
        let fenestrations =
            SolarSurface::make_fenestrations(&model.fenestrations, state, n_solar_rays)?;

        // build scene
        let mut solar_scene = Scene::from_simple_model(model, Wavelengths::Solar)?;
        solar_scene.build_accelerator();

        // calculator
        let solar_dc_factory = DCFactory {
            max_depth: 0,
            n_ambient_samples: *options.solar_ambient_divitions().unwrap(),
            reinhart: ReinhartSky::new(mf),
            ..DCFactory::default()
        };

        // calculate
        let front_surfaces_dc =
            SolarSurface::calc_solar_dc_matrix(&surfaces, &solar_scene, &solar_dc_factory, true);

        let back_surfaces_dc =
            SolarSurface::calc_solar_dc_matrix(&surfaces, &solar_scene, &solar_dc_factory, false);

        let front_fenestrations_dc = SolarSurface::calc_solar_dc_matrix(
            &fenestrations,
            &solar_scene,
            &solar_dc_factory,
            true,
        );

        let back_fenestrations_dc = SolarSurface::calc_solar_dc_matrix(
            &fenestrations,
            &solar_scene,
            &solar_dc_factory,
            false,
        );

        let front_surfaces_view_factors: Vec<IRViewFactorSet> = surfaces
            .iter()
            .map(|s| s.calc_view_factors(&solar_scene, true))
            .collect();
        let back_surfaces_view_factors: Vec<IRViewFactorSet> = surfaces
            .iter()
            .map(|s| s.calc_view_factors(&solar_scene, false))
            .collect();
        let front_fenestrations_view_factors: Vec<IRViewFactorSet> = fenestrations
            .iter()
            .map(|s| s.calc_view_factors(&solar_scene, true))
            .collect();
        let back_fenestrations_view_factors: Vec<IRViewFactorSet> = fenestrations
            .iter()
            .map(|s| s.calc_view_factors(&solar_scene, false))
            .collect();

        Ok(Self {
            front_surfaces_view_factors,
            back_surfaces_view_factors,
            front_fenestrations_view_factors,
            back_fenestrations_view_factors,
            front_surfaces_dc,
            back_surfaces_dc,
            front_fenestrations_dc,
            back_fenestrations_dc,
        })
    }
}
