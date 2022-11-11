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
use crate::Float;
use calendar::Date;
use communication_protocols::{ErrorHandling, MetaOptions, SimulationModel};
use matrix::Matrix;
use simple_model::{Boundary, SimpleModel, SimulationState, SimulationStateHeader, SolarOptions};
use solar::{PerezSky, SkyUnits, Solar};
use std::borrow::Borrow;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use weather::{CurrentWeather, Weather};

use crate::optical_info::OpticalInfo;

/// The main model
pub struct SolarModel {
    // /// The scene that makes up this model from a lighting point of view.
    // lighting_scene: Scene,

    // Workplanes
    /// The scene that makes up this model from a radiation point of view.
    // solar_scene: Scene,

    // surfaces: Vec<SolarSurface>,

    /// The optical information from the model, containing
    /// DC matrices and view factors
    optical_info: OpticalInfo,

    /// The calculator for solar position and other solar variables
    solar: Solar,

    /// The MF discretization scheme for the sky.
    solar_sky_discretization: usize,
}

impl SolarModel {
    /// This function makes the IR heat transfer Zero... we will try to fix this soon enough,
    /// just not now    
    fn update_ir_radiation(
        &self,
        weather_data: &CurrentWeather,
        model: &SimpleModel,
        state: &mut SimulationState,
    ) -> Result<(), String> {
        pub const SIGMA: crate::Float = 5.670374419e-8;

        fn ir(temp: Float, emissivity: Float) -> Float {
            emissivity * SIGMA * (temp + 273.15).powi(4)
        }

        let db = match weather_data.dry_bulb_temperature {
            Some(v) => v,
            None => return Err("Cannot calculate IR radiation without Dry Bulb temperature".into()),
        };
        let horizontal_ir = match weather_data.horizontal_infrared_radiation_intensity {
            Some(v) => v,
            None => weather_data.derive_horizontal_ir()?,
        };

        let iter = model.surfaces.iter().enumerate();

        for (index, surface) in iter {
            // Deal with front
            if let Ok(b) = surface.front_boundary() {
                if let Boundary::Space { .. } = b {
                    // let temp = space.dry_bulb_temperature(state).unwrap_or_else(|| 22.);
                    let temp = surface.first_node_temperature(state).unwrap_or( 22.);
                    surface.set_front_ir_irradiance(state, ir(temp, 1.0));
                } // else is ground... ignore
            } else {
                // outdoor
                let view_factors = &self.optical_info.front_surfaces_view_factors[index];
                let ground_other = (view_factors.ground + view_factors.air) * ir(db, 1.0);
                let sky = view_factors.sky * horizontal_ir;
                surface.set_front_ir_irradiance(state, ground_other + sky);
            }

            // Deal with Back
            if let Ok(b) = surface.back_boundary() {
                if let Boundary::Space { .. } = b {
                    // let temp = space.dry_bulb_temperature(state).unwrap_or_else(|| 22.);
                    let temp = surface.last_node_temperature(state).unwrap_or(22.);
                    surface.set_back_ir_irradiance(state, ir(temp, 1.0));
                } // else is ground... ignore
            } else {
                // outdoor
                let view_factors = &self.optical_info.back_surfaces_view_factors[index];
                let ground_other = (view_factors.ground + view_factors.air) * ir(db, 1.0);
                let sky = view_factors.sky * horizontal_ir;
                surface.set_back_ir_irradiance(state, ground_other + sky);
            }
        }

        let iter = model.fenestrations.iter().enumerate();
        for (index, surface) in iter {
            // Deal with front
            if let Ok(b) = surface.front_boundary() {
                if let Boundary::Space { .. } = b {
                    // let temp = space.dry_bulb_temperature(state).unwrap_or_else(|| 22.);
                    let temp = surface.first_node_temperature(state).unwrap_or( 22.);
                    surface.set_front_ir_irradiance(state, ir(temp, 1.0));
                } // else is ground... ignore
            } else {
                // outdoor
                let view_factors = &self.optical_info.front_fenestrations_view_factors[index];
                let ground_other = (view_factors.ground + view_factors.air) * ir(db, 1.0);
                let sky = view_factors.sky * horizontal_ir;
                surface.set_front_ir_irradiance(state, ground_other + sky);
            }

            // Deal with Back
            if let Ok(b) = surface.back_boundary() {
                if let Boundary::Space { .. } = b {
                    // let temp = space.dry_bulb_temperature(state).unwrap_or_else(|| 22.);
                    let temp = surface.last_node_temperature(state).unwrap_or( 22.);
                    surface.set_back_ir_irradiance(state, ir(temp, 1.0));
                } // else is ground... ignore
            } else {
                // outdoor
                let view_factors = &self.optical_info.back_fenestrations_view_factors[index];
                let ground_other = (view_factors.ground + view_factors.air) * ir(db, 1.0);
                let sky = view_factors.sky * horizontal_ir;
                surface.set_back_ir_irradiance(state, ground_other + sky);
            }
        }

        Ok(())
    }

    fn update_solar_radiation(
        &self,
        date: Date,
        weather_data: CurrentWeather,
        model: &SimpleModel,
        state: &mut SimulationState,
    ) {
        let direct_normal_irrad = weather_data
            .direct_normal_radiation
            .expect("Missing data for direct normal irradiance");
        let diffuse_horizontal_irrad = weather_data
            .diffuse_horizontal_radiation
            .expect("Missing data for diffuse horizontal");

        let is_day = direct_normal_irrad + diffuse_horizontal_irrad >= 1e-4;
        let vec = if is_day {
            // Build sky vector
            let albedo = 0.2;
            let add_sky = true;
            let add_sun = true;
            let units = SkyUnits::Solar;
            PerezSky::gen_sky_vec(
                self.solar_sky_discretization,
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
        if !self.optical_info.front_surfaces_dc.is_empty() {
            if is_day {
                let solar_irradiance = &self.optical_info.front_surfaces_dc * &vec;
                for (i, s) in model.surfaces.iter().enumerate() {
                    // Average of the period
                    let mut v = solar_irradiance.get(i, 0).unwrap();
                    if v < 0.0 {
                        v = 0.0
                    }
                    let old_v = s.front_incident_solar_irradiance(state).unwrap();
                    s.set_front_incident_solar_irradiance(state, (v + old_v) / 2.);
                }
            } else {
                for s in model.surfaces.iter() {
                    s.set_front_incident_solar_irradiance(state, 0.0);
                }
            }
        }
        if !self.optical_info.back_surfaces_dc.is_empty() {
            if is_day {
                let solar_irradiance = &self.optical_info.back_surfaces_dc * &vec;
                for (i, s) in model.surfaces.iter().enumerate() {
                    // Average of the period
                    let mut v = solar_irradiance.get(i, 0).unwrap();
                    if v < 0.0 {
                        v = 0.0
                    }
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
        if !self.optical_info.front_fenestrations_dc.is_empty() {
            if is_day {
                let solar_irradiance = &self.optical_info.front_fenestrations_dc * &vec;
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
        if !self.optical_info.back_fenestrations_dc.is_empty() {
            if is_day {
                let solar_irradiance = &self.optical_info.back_fenestrations_dc * &vec;
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
    }
}

impl ErrorHandling for SolarModel {
    fn module_name() -> &'static str {
        "Solar Model"
    }
}

impl SimulationModel for SolarModel {
    type Type = Self;
    type OptionType = SolarOptions;
    fn new<M: Borrow<SimpleModel>>(
        meta_options: &MetaOptions,
        options: SolarOptions,
        model: M,
        state: &mut SimulationStateHeader,
        _n: usize,
    ) -> Result<Self::Type, String> {
        let model = model.borrow();
        // Make OpticalInfo, or read, as needed
        let optical_info = if let Ok(path_str) = options.optical_data_path() {
            let path = Path::new(path_str);
            if path.exists() {
                // read from file
                assert!(
                    path.is_file(),
                    "Path '{}' is not a file",
                    path.to_str().unwrap()
                );
                let data = match std::fs::read_to_string(path) {
                    Ok(v) => v,
                    Err(_) => {
                        return Err(format!("Unable to read optical_info file '{}'", path_str))
                    }
                };
                let info: OpticalInfo = match serde_json::from_str(&data) {
                    Ok(v) => v,
                    Err(_) => {
                        return Err(format!(
                            "Unable to patse optical_info object in file '{}'",
                            path_str
                        ))
                    }
                };

                info
            } else {
                // write into file
                let info = OpticalInfo::new(&options, model, state)?;
                let s = serde_json::to_value(&info).unwrap();
                let mut file = File::create(path).unwrap();
                writeln!(&mut file, "{}", s).unwrap();
                info
            }
        } else {
            // Forced calculation... not store
            OpticalInfo::new(&options, model, state)?
        };

        // Create the Solar object
        let latitude = meta_options.latitude;
        let longitude = -meta_options.longitude;
        let standard_meridian = -meta_options.standard_meridian;
        let solar = Solar::new(latitude, longitude, standard_meridian);

        // derive MF
        let (.., ncols) = optical_info.back_surfaces_dc.size();
        if ncols == 0 {
            return Err(
                "optical data is corrupt: daylight coefficient matrix has zero columns.".to_string()
            );
        }
        let mut mf = 1;
        loop {
            if mf >= 9 || ncols == 0 {
                return Err(format!("sky discretization seems to be too high ({mf}... If this is a bug, please report it!"));
            }
            if solar::ReinhartSky::n_bins(mf) == ncols {
                break;
            } else {
                mf += 1;
            }
        }

        Ok(Self {
            optical_info,
            solar,
            solar_sky_discretization: mf,
        })
    }

    fn march<W: Weather, M: Borrow<SimpleModel>>(
        &self,
        date: Date,
        weather: &W,
        model: M,
        state: &mut SimulationState,
    ) -> Result<(), String> {
        let model = model.borrow();
        // Handle the solar part

        let weather_data = weather.get_weather_data(date);

        self.update_ir_radiation(&weather_data, model, state)?;
        self.update_solar_radiation(date, weather_data, model, state);

        // return
        Ok(())
        // unimplemented!()
    }
}

#[cfg(test)]
mod testing {
    use super::*;

    #[test]
    fn test_model_mf() {
        // cleanup
        let optical_data_path = "./tests/wall/optical_data.json";
        let path = Path::new(optical_data_path);
        if path.exists() {
            std::fs::remove_file(path).unwrap();
        }

        // Step 1: create the model, write the optical info data.
        let meta_options = MetaOptions {
            latitude: -33.,
            longitude: 72.,
            standard_meridian: 70.,
            elevation: 0.0,
        };
        let (model, mut state_header) = SimpleModel::from_file("./tests/wall/wall.spl").unwrap();
        let mut solar_options = model.solar_options.clone().unwrap();
        solar_options.set_optical_data_path(optical_data_path.to_string());

        let light_model =
            SolarModel::new(&meta_options, solar_options, &model, &mut state_header, 4).unwrap();
        assert_eq!(light_model.solar_sky_discretization, 1); //this comes in the model

        // Step 2: Run it again, with a different option
        let meta_options = MetaOptions {
            latitude: -33.,
            longitude: 72.,
            standard_meridian: 70.,
            elevation: 0.0,
        };
        let mut solar_options = model.solar_options.clone().unwrap();
        solar_options.set_solar_sky_discretization(2);
        solar_options.set_optical_data_path(optical_data_path.to_string());
        let light_model =
            SolarModel::new(&meta_options, solar_options, &model, &mut state_header, 4).unwrap();
        assert_eq!(light_model.solar_sky_discretization, 1); //this comes from the optical data.

        // cleanup
        let path = Path::new(optical_data_path);
        if path.exists() {
            std::fs::remove_file(path).unwrap();
        }
    }
}
