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
use calendar::Date;
use communication_protocols::{ErrorHandling, MetaOptions, SimulationModel};
use matrix::Matrix;
use rendering::{DCFactory, Scene, Wavelengths};
use simple_model::{SimpleModel, SimulationState, SimulationStateHeader, SolarOptions};
use solar::ReinhartSky;
use solar::{PerezSky, SkyUnits, Solar};
use weather::Weather;



/// The main model
pub struct SolarModel {
    // /// The scene that makes up this model from a lighting point of view.
    // lighting_scene: Scene,

    // Workplanes
    /// The scene that makes up this model from a radiation point of view.
    // solar_scene: Scene,

    // surfaces: Vec<SolarSurface>,

    /// The Daylight Coefficients matrix for the front-side of the  surfaces in the scene
    front_surfaces_dc: Matrix,

    /// The Daylight Coefficients matrix for the back-side of the  surfaces in the scene
    back_surfaces_dc: Matrix,

    // fenestrations: Vec<SolarSurface>,
    /// The Daylight Coefficients matrix for the front-side of the  fenestrations in the scene
    front_fenestrations_dc: Matrix,

    /// The Daylight Coefficients matrix for the back-side of the fenestrations in the scene
    back_fenestrations_dc: Matrix,

    /// The options for the model.
    options: SolarOptions,

    /// The calculator for solar position and other solar variables
    solar: Solar,
}

impl ErrorHandling for SolarModel {
    fn module_name() -> &'static str {
        "Solar Model"
    }
}

impl SimulationModel for SolarModel {
    type Type = Self;
    type OptionType = SolarOptions;
    fn new(
        meta_options: &MetaOptions,
        options: SolarOptions,
        model: &SimpleModel,
        state: &mut SimulationStateHeader,
        _n: usize,
    ) -> Result<Self::Type, String> {
        let latitude = meta_options.latitude;
        let longitude = -meta_options.longitude;
        let standard_meridian = -meta_options.standard_meridian;
        let solar = Solar::new(latitude, longitude, standard_meridian);

        // let lighting_scene = Scene::from_simple_model(&model, Wavelengths::Visible);
        // lighting_scene.build_accelerator();

        /* *********************** */
        /* PROCESS SOLAR RADIATION */
        /* *********************** */
        let mut solar_scene = Scene::from_simple_model(model, Wavelengths::Solar);
        solar_scene.build_accelerator();
        let mf = *options.solar_sky_discretization().unwrap();
        let n_solar_rays = *options.n_solar_irradiance_points().unwrap();

        let solar_dc_factory = DCFactory {
            max_depth: 0,
            n_ambient_samples: *options.solar_ambient_divitions().unwrap(),
            reinhart: ReinhartSky::new(mf),
            ..DCFactory::default()
        };

        // Create Surfaces
        let surfaces = SolarSurface::make_surfaces(&model.surfaces, state, n_solar_rays);
        let path = match options.front_surfaces_solar_irradiance_matrix() {
            Ok(e) => Some(e),
            Err(_e) => None,
        };
        let front_surfaces_dc = SolarSurface::get_front_solar_dc_matrix(
            &surfaces,
            path,
            &solar_scene,
            &solar_dc_factory,
        )?;
        let path = match options.back_surfaces_solar_irradiance_matrix() {
            Ok(e) => Some(e),
            Err(_e) => None,
        };
        let back_surfaces_dc = SolarSurface::get_back_solar_dc_matrix(
            &surfaces,
            path,
            &solar_scene,
            &solar_dc_factory,
        )?;

        // Process Fenestrations
        let path = match options.front_fenestrations_solar_irradiance_matrix() {
            Ok(e) => Some(e),
            Err(_e) => None,
        };
        let fenestrations =
            SolarSurface::make_fenestrations(&model.fenestrations, state, n_solar_rays);
        let front_fenestrations_dc = SolarSurface::get_front_solar_dc_matrix(
            &fenestrations,
            path,
            &solar_scene,
            &solar_dc_factory,
        )?;

        let path = match options.back_fenestrations_solar_irradiance_matrix() {
            Ok(e) => Some(e),
            Err(_e) => None,
        };
        let back_fenestrations_dc = SolarSurface::get_back_solar_dc_matrix(
            &fenestrations,
            path,
            &solar_scene,
            &solar_dc_factory,
        )?;

        Ok(Self {
            options,
            // solar_scene,
            // surfaces,
            front_surfaces_dc,
            back_surfaces_dc,
            // fenestrations,
            front_fenestrations_dc,
            back_fenestrations_dc,
            solar,
        })
    }

    fn march(
        &self,
        date: Date,
        weather: &dyn Weather,
        model: &SimpleModel,
        state: &mut SimulationState,
    ) -> Result<(), String> {
        // Handle the solar part

        let weather_data = weather.get_weather_data(date);

        // let date = Date::from_day_of_year(self.solar.unwrap_solar_time(Time::Solar(date.day_of_year())));

        let direct_normal_irrad = weather_data.direct_normal_radiation.unwrap();
        let diffuse_horizontal_irrad = weather_data.diffuse_horizontal_radiation.unwrap();
        let is_day = direct_normal_irrad + diffuse_horizontal_irrad >= 1e-4;
        let vec = if is_day {
            // Build sky vector
            let albedo = 0.2;
            let add_sky = true;
            let add_sun = true;
            let units = SkyUnits::Solar;
            PerezSky::gen_sky_vec(
                *self.options.solar_sky_discretization().unwrap(),
                &self.solar,
                date,
                weather_data,
                units,
                albedo,
                add_sky,
                add_sun,
            )
            .unwrap()
        } else {
            Matrix::empty()
        };

        // Process Solar Irradiance in Surfaces
        if !self.front_surfaces_dc.is_empty() {
            if is_day {
                let solar_irradiance = &self.front_surfaces_dc * &vec;
                if solar_irradiance.get(0, 0).unwrap() < 0.0 {
                    dbg!(solar_irradiance.get(0, 0).unwrap());
                }
                for (i, s) in model.surfaces.iter().enumerate() {
                    // Average of the period
                    let v = solar_irradiance.get(i, 0).unwrap();
                    let old_v = s.front_incident_solar_irradiance(state).unwrap();
                    s.set_front_incident_solar_irradiance(state, (v + old_v) / 2.);
                }
            } else {
                for s in model.surfaces.iter() {
                    s.set_front_incident_solar_irradiance(state, 0.0);
                }
            }
        }
        if !self.back_surfaces_dc.is_empty() {
            if is_day {
                let solar_irradiance = &self.back_surfaces_dc * &vec;
                for (i, s) in model.surfaces.iter().enumerate() {
                    // Average of the period
                    let v = solar_irradiance.get(i, 0).unwrap();
                    let old_v = s.back_incident_solar_irradiance(state).unwrap();
                    s.set_back_incident_solar_irradiance(state, (v + old_v) / 2.);
                }
            } else {
                for s in model.surfaces.iter() {
                    s.set_front_incident_solar_irradiance(state, 0.0);
                }
            }
        }

        // Process Solar Irradiance in Fenestration
        if !self.front_fenestrations_dc.is_empty() {
            if is_day {
                let solar_irradiance = &self.front_fenestrations_dc * &vec;
                for (i, s) in model.fenestrations.iter().enumerate() {
                    // Average of the period
                    let v = solar_irradiance.get(i, 0).unwrap();
                    let old_v = s.front_incident_solar_irradiance(state).unwrap();
                    s.set_front_incident_solar_irradiance(state, (v + old_v) / 2.);
                }
            } else {
                for s in model.fenestrations.iter() {
                    s.set_front_incident_solar_irradiance(state, 0.0);
                }
            }
        }
        if !self.back_fenestrations_dc.is_empty() {
            if is_day {
                let solar_irradiance = &self.back_fenestrations_dc * &vec;
                for (i, s) in model.fenestrations.iter().enumerate() {
                    // Average of the period
                    let v = solar_irradiance.get(i, 0).unwrap();
                    let old_v = s.back_incident_solar_irradiance(state).unwrap();
                    s.set_back_incident_solar_irradiance(state, (v + old_v) / 2.);
                }
            } else {
                for s in model.fenestrations.iter() {
                    s.set_front_incident_solar_irradiance(state, 0.0);
                }
            }
        }

        // return
        Ok(())
        // unimplemented!()
    }
}
