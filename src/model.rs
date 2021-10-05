use weather::Weather;
use simple_model::model::SimpleModel;
use calendar::date::Date;
use communication_protocols::error_handling::ErrorHandling;
use communication_protocols::simulation_model::SimulationModel;
use simple_model::simulation_state::{SimulationState,SimulationStateHeader};
use simple_model::simulation_state_element::SimulationStateElement;

pub struct SolarModel ();

impl ErrorHandling for SolarModel {
    fn module_name() -> &'static str {
        "Placeholder Solar Model"
    }
}

impl SimulationModel for SolarModel {
    type Type = Self;

    fn new(model : &SimpleModel, state: &mut SimulationStateHeader, _n: usize)->Result<Self::Type,String>{

        for (i,space) in model.spaces.iter().enumerate(){
            // Initialize at night.... this will change right away because
            // light is quasi-static (does not depend on the apast)
            let index = state.push( SimulationStateElement::SpaceBrightness(i), 0.0);
            space.set_brightness_index(index);
        }

        // We could do the same thing with SolarRadiation over walls, but the 
        // current THermalModule does not support that.

        Ok(Self())
    }

    fn march(&self, date: Date, weather: &dyn Weather, model: &SimpleModel, state: &mut SimulationState)->Result<(),String>{
        let current_weather = weather.get_weather_data(date);

        let direct_normal_radiation = current_weather.direct_normal_radiation.unwrap();
        let global_horizontal_radiation = current_weather.global_horizontal_radiation.unwrap();

        for space in model.spaces.iter(){
            
            let brightness = 0.4*global_horizontal_radiation + 0.5 * direct_normal_radiation;
            space.set_brightness(state, brightness)
                  
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
