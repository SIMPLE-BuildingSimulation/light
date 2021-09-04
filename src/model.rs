use weather::Weather;
use building_model::building::Building;
use calendar::date::Date;
use communication_protocols::error_handling::ErrorHandling;
use communication_protocols::simulation_model::SimulationModel;
use building_model::simulation_state::SimulationState;
use building_model::simulation_state_element::SimulationStateElement;

pub struct SolarModel ();

impl ErrorHandling for SolarModel {
    fn module_name() -> &'static str {
        "Placeholder Solar Model"
    }
}

impl SimulationModel for SolarModel {
    type Type = Self;

    fn new(building : &Building, state: &mut SimulationState, _n: usize)->Result<Self::Type,String>{

        for space in &building.spaces{
            // Initialize at night.... this will change right away because
            // light is quasi-static (does not depend on the apast)
            state.push( SimulationStateElement::SpaceBrightness(space.index().unwrap(), 0.0));
        }

        // We could do the same thing with SolarRadiation over walls, but the 
        // current THermalModule does not support that.

        Ok(Self())
    }

    fn march(&self, date: Date, weather: &dyn Weather, building: &Building, state: &mut SimulationState)->Result<(),String>{
        let current_weather = weather.get_weather_data(date);

        let direct_normal_radiation = current_weather.direct_normal_radiation.unwrap();
        let global_horizontal_radiation = current_weather.global_horizontal_radiation.unwrap();

        for (space_index,space) in building.spaces.iter().enumerate(){
            let i = space.brightness_index().unwrap();

            state.update_value(i, 
                SimulationStateElement::SpaceBrightness(space_index, 0.4*global_horizontal_radiation + 0.5 * direct_normal_radiation)
            )      
        }

        Ok(())
    }


}


/***********/
/* TESTING */
/***********/



#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
