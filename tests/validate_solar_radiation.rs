use calendar::Date;
use communication_protocols::{MetaOptions, SimulationModel};
use schedule::ScheduleConstant;
use simple_model::SolarOptions;
use simple_test_models::*;
use light::{Float, SolarModel};
use validate::{valid, SeriesValidator, Validate, Validator};
use weather::SyntheticWeather;

fn get_validator(expected: Vec<f64>, found: Vec<f64>) -> Box<SeriesValidator> {
    Box::new(SeriesValidator {
        x_label: Some("Timestep".into()),
        y_label: Some("Solar Irradiance".into()),
        y_units: Some("W/m2"),
        expected_legend: Some("EnergyPlus"),
        found_legend: Some("SIMPLE"),
        expected,
        found,
        ..validate::SeriesValidator::default()
    })
}

fn get_expected(city: &str, orientation: &str) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let path = format!("./tests/{city}_{orientation}/eplusout.csv");
    let cols = validate::from_csv(&path, &[2, 3, 4]);
    let diffuse_horizontal_rad = &cols[0];
    let direct_normal_rad = &cols[1]; 
    let incident_solar_radiation = &cols[2]; 

    (
        incident_solar_radiation.clone(),
        diffuse_horizontal_rad.clone(),
        direct_normal_rad.clone(),
    )
}

fn get_simple_results(
    city: &str,
    diffuse_horizontal_rad: Vec<f64>,
    direct_normal_rad: Vec<f64>,
    orientation: f64,
) -> Vec<f64> {
    let (lat, lon, std_mer): (Float, Float, Float) = match city.as_bytes() {
        b"wellington" => (-41.3, 174.78, 180.),
        b"barcelona" => (41.28, 2.07, 15.), // ??? GMT + 1
        _ => panic!("Unsupported city '{}'", city),
    };

    let meta_options = MetaOptions {
        latitude: lat.to_radians(),
        longitude: lon.to_radians(),
        standard_meridian: std_mer.to_radians(),
        elevation: 0.0,
    };

    
    let zone_volume = 600.;

    let (simple_model, mut state_header) = get_single_zone_test_building(        
        &SingleZoneTestBuildingOptions {
            zone_volume,
            surface_width: 20.,
            surface_height: 3.,
            construction: vec![TestMat::Concrete(0.2)],
            orientation,
            ..Default::default()
        },
    );

    // Finished model the SimpleModel
    let mut options = SolarOptions::new();
    options
        .set_n_solar_irradiance_points(10)
        .set_solar_ambient_divitions(3000)
        .set_solar_sky_discretization(1);

    let n: usize = 20;
    let solar_model =
        SolarModel::new(&meta_options, options, &simple_model, &mut state_header, n).unwrap();
    let mut state = state_header.take_values().unwrap();
    let mut date = Date {
        month: 1,
        day: 1,
        hour: 0.5,
    };
    let mut ret = Vec::with_capacity(diffuse_horizontal_rad.len());
    for (diffuse_horizontal, direct_normal) in
        diffuse_horizontal_rad.iter().zip(direct_normal_rad.iter())
    {
        // Set outdoor temp
        let mut weather = SyntheticWeather::default();
        weather.direct_normal_radiation = Box::new(ScheduleConstant::new(*direct_normal));
        weather.diffuse_horizontal_radiation = Box::new(ScheduleConstant::new(*diffuse_horizontal));
        weather.dew_point_temperature = Box::new(ScheduleConstant::new(11.)); // 11C is what Radiance uses by default.
        weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(21.)); // should be irrelevant
        weather.opaque_sky_cover = Box::new(ScheduleConstant::new(0.)); // should be irrelevant

        let surface = &simple_model.surfaces[0];

        // March
        solar_model
            .march(date, &weather, &simple_model, &mut state)
            .unwrap();

        let front_radiation = surface.front_incident_solar_irradiance(&state).unwrap();
        ret.push(front_radiation);
        // let back_radiation = surface.back_incident_solar_irradiance(&state).unwrap();

        // let lat = -33.38 * PI / 180.;
        // let lon = 70.78 * PI / 180.;
        // let std_mer = 60. * PI / 180.;

        // let solar = Solar::new(lat, lon, std_mer);
        // let day = Time::Standard(date.day_of_year());
        // let err = (front_radiation - exp_radiation).abs();

        // println!("{},{},{}", front_radiation, exp_radiation, err);

        // assert!( (found_radiation-exp_radiation).abs() < 0.4, "found_temp = {}, exp_temp = {} ,  error = {}", found_radiation, exp_radiation, (found_radiation - *exp_radiation).abs() );

        // Advance
        date.add_hours(1. / n as f64);
        // assert!(false)
    }
    ret
}

fn barcelona(validator: &mut Validator) {
    const CITY: &'static str = "barcelona";

    #[valid(Exterior Incident Solar Radiation - Barcelona, South)]
    fn validate_barcelona_south() -> Box<dyn Validate> {
        let (expected, diffuse_horizontal_rad, direct_normal_rad) = get_expected(CITY, "south");
        let found = get_simple_results(CITY, diffuse_horizontal_rad, direct_normal_rad, 0.0);
        get_validator(expected, found)
    }

    #[valid(Exterior Incident Solar Radiation - Barcelona, North)]
    fn validate_barcelona_north() -> Box<dyn Validate> {
        let (expected, diffuse_horizontal_rad, direct_normal_rad) = get_expected(CITY, "north");
        let found = get_simple_results(CITY, diffuse_horizontal_rad, direct_normal_rad, 180.0);
        get_validator(expected, found)
    }

    #[valid(Exterior Incident Solar Radiation - Barcelona, West)]
    fn validate_barcelona_west() -> Box<dyn Validate> {
        let (expected, diffuse_horizontal_rad, direct_normal_rad) = get_expected(CITY, "west");
        let found = get_simple_results(CITY, diffuse_horizontal_rad, direct_normal_rad, 90.0);
        get_validator(expected, found)
    }

    #[valid(Exterior Incident Solar Radiation - Barcelona, East)]
    fn validate_barcelona_east() -> Box<dyn Validate> {
        let (expected, diffuse_horizontal_rad, direct_normal_rad) = get_expected(CITY, "east");
        let found = get_simple_results(CITY, diffuse_horizontal_rad, direct_normal_rad, -90.);
        get_validator(expected, found)
    }

    validator.push(validate_barcelona_south());
    validator.push(validate_barcelona_north());
    validator.push(validate_barcelona_west());
    validator.push(validate_barcelona_east());
}

fn wellington(validator: &mut Validator) {
    const CITY: &'static str = "wellington";

    #[valid(Exterior Incident Solar Radiation - Wellington, South)]
    fn validate_wellington_south() -> Box<dyn Validate> {
        let (expected, diffuse_horizontal_rad, direct_normal_rad) = get_expected(CITY, "south");
        let found = get_simple_results(CITY, diffuse_horizontal_rad, direct_normal_rad, 0.0);
        get_validator(expected, found)
    }

    #[valid(Exterior Incident Solar Radiation - Wellington, North)]
    fn validate_wellington_north() -> Box<dyn Validate> {
        let (expected, diffuse_horizontal_rad, direct_normal_rad) = get_expected(CITY, "north");
        let found = get_simple_results(CITY, diffuse_horizontal_rad, direct_normal_rad, 180.0);
        get_validator(expected, found)
    }

    #[valid(Exterior Incident Solar Radiation - Wellington, West)]
    fn validate_wellington_west() -> Box<dyn Validate> {
        let (expected, diffuse_horizontal_rad, direct_normal_rad) = get_expected(CITY, "west");
        let found = get_simple_results(CITY, diffuse_horizontal_rad, direct_normal_rad, 90.0);
        get_validator(expected, found)
    }

    #[valid(Exterior Incident Solar Radiation - Wellington, East)]
    fn validate_wellington_east() -> Box<dyn Validate> {
        let (expected, diffuse_horizontal_rad, direct_normal_rad) = get_expected(CITY, "east");
        let found = get_simple_results(CITY, diffuse_horizontal_rad, direct_normal_rad, -90.);
        get_validator(expected, found)
    }

    validator.push(validate_wellington_south());
    validator.push(validate_wellington_north());
    validator.push(validate_wellington_west());
    validator.push(validate_wellington_east());
}

#[test]
fn validate_solar_radiation() {
    // cargo test --package light --test validate_solar_radiation -- validate_solar_radiation --exact --nocapture
    let mut validator = Validator::new(
        "Validate Solar Radiation",
        "./docs/validation/incident_solar_radiation.html",
    );

    barcelona(&mut validator);
    wellington(&mut validator);

    validator.validate().unwrap();
}
