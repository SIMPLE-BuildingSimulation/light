var searchIndex = JSON.parse('{\
"light":{"doc":"This is SIMPLE’s Light and Solar calculation module. It …","t":[6,3,3,17,2,12,12,12,12,12,12,12,12,12,12,0,0,12,0,3,11,11,11,11,11,11,11,12,12,12,11,11,11,11,11,11,3,3,12,12,12,12,12,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,11,11,11,11,11,12,11,11,11,11,11,11,11,11,11,11,18,3,11,11,11,11,11,11,11,5,11,11,11,11,12,12,11,11,11,11,11],"n":["Float","IRViewFactorSet","OpticalInfo","PI","SolarModel","air","back_fenestrations_dc","back_fenestrations_view_factors","back_surfaces_dc","back_surfaces_view_factors","front_fenestrations_dc","front_fenestrations_view_factors","front_surfaces_dc","front_surfaces_view_factors","ground","model","optical_info","sky","solar_surface","SolarModel","borrow","borrow_mut","from","into","march","module_name","new","optical_info","solar","solar_sky_discretization","try_from","try_into","type_id","update_ir_radiation","update_solar_radiation","vzip","IRViewFactorSet","OpticalInfo","air","back_fenestrations_dc","back_fenestrations_view_factors","back_surfaces_dc","back_surfaces_view_factors","borrow","borrow","borrow_mut","borrow_mut","clone","clone","clone_into","clone_into","default","deserialize","deserialize","fmt","from","from","front_fenestrations_dc","front_fenestrations_view_factors","front_surfaces_dc","front_surfaces_view_factors","ground","into","into","new","serialize","serialize","sky","to_owned","to_owned","try_from","try_from","try_into","try_into","type_id","type_id","vzip","vzip","DELTA","SolarSurface","back_rays","borrow","borrow_mut","calc_solar_dc_matrix","calc_view_factors","from","front_rays","get_sampler","into","make_fenestrations","make_surfaces","new","normal","points","solar_irradiance","try_from","try_into","type_id","vzip"],"q":["light","","","","","","","","","","","","","","","","","","","light::model","","","","","","","","","","","","","","","","","light::optical_info","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","light::solar_surface","","","","","","","","","","","","","","","","","","","",""],"d":["The kind of Floating point number used in the library… …","A set of view factors as seen by a <code>ThermalSurface</code>.","Information about the solar radiation and other optical …","Well, Pi.","","The fraction of the view that corresponds to other objects …","The Daylight Coefficients matrix for the back-side of the …","The <code>IRViewFactorSet</code> for the back side of each fenestration","The Daylight Coefficients matrix for the back-side of the  …","The <code>IRViewFactorSet</code> for the back side of each surface","The Daylight Coefficients matrix for the front-side of the …","The <code>IRViewFactorSet</code> for the front side of each fenestration","The Daylight Coefficients matrix for the front-side of the …","The <code>IRViewFactorSet</code> for the front side of each surface","The fraction of the view that corresponds to the ground","The main export of this module: A Simulation Model for …","","The fraction of the view that corresponds to the sky","","The main model","","","Returns the argument unchanged.","Calls <code>U::from(self)</code>.","","","","The scene that makes up this model from a radiation point …","The calculator for solar position and other solar variables","The MF discretization scheme for the sky.","","","","This function makes the IR heat transfer Zero… we will …","","","A set of view factors as seen by a <code>ThermalSurface</code>.","Information about the solar radiation and other optical …","The fraction of the view that corresponds to other objects …","The Daylight Coefficients matrix for the back-side of the …","The <code>IRViewFactorSet</code> for the back side of each fenestration","The Daylight Coefficients matrix for the back-side of the  …","The <code>IRViewFactorSet</code> for the back side of each surface","","","","","","","","","","","","","Returns the argument unchanged.","Returns the argument unchanged.","The Daylight Coefficients matrix for the front-side of the …","The <code>IRViewFactorSet</code> for the front side of each fenestration","The Daylight Coefficients matrix for the front-side of the …","The <code>IRViewFactorSet</code> for the front side of each surface","The fraction of the view that corresponds to the ground","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Calculates the new OpticalInformation ","","","The fraction of the view that corresponds to the sky","","","","","","","","","","","Offset for the starting point of the rays.","Structure that can help calculate solar radiation","Gets the back rays of a surface","","","Receives a list of <code>SolarSurface</code> objects as well as the …","Calculates an <code>IRViewFactorSet</code> for this surface","Returns the argument unchanged.","Gets the front rays of a surface","","Calls <code>U::from(self)</code>.","Builds a set of SolarSurfaces from Fenestrations","Builds a set of SolarSurfaces from Surfaces","","","","Calculates the Daylight Coefficient matrix for the front …","","","",""],"i":[0,0,0,0,0,14,15,15,15,15,15,15,15,15,14,0,0,14,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,14,15,15,15,15,14,15,14,15,14,15,14,15,14,14,15,14,14,15,15,15,15,15,14,14,15,15,14,15,14,14,15,14,15,14,15,14,15,14,15,18,0,18,18,18,18,18,18,18,0,18,18,18,18,18,18,18,18,18,18,18],"f":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,[[]],[[]],[[]],[[]],[[1,2,3,4],[[6,[5]]]],[[],7],[[8,9,3,10,11],[[6,[5]]]],0,0,0,[[],6],[[],6],[[],12],[[1,13,3,4],[[6,[5]]]],[[1,2,13,3,4]],[[]],0,0,0,0,0,0,0,[[]],[[]],[[]],[[]],[14,14],[15,15],[[]],[[]],[[],14],[[],[[6,[14]]]],[[],[[6,[15]]]],[[14,16],17],[[]],[[]],0,0,0,0,0,[[]],[[]],[[9,3,10],15],[14,6],[15,6],0,[[]],[[]],[[],6],[[],6],[[],6],[[],6],[[],12],[[],12],[[]],[[]],0,0,[18,[[20,[19]]]],[[]],[[]],[[21,22,23],24],[[18,21,23],14],[[]],[18,[[20,[19]]]],[[[20,[25]]],26],[[]],[[10,11],[[20,[18]]]],[[10,11],[[20,[18]]]],[[11,27],18],0,0,[[18,21,22],24],[[],6],[[],6],[[],12],[[]]],"p":[[3,"SolarModel"],[3,"Date"],[3,"SimpleModel"],[6,"SimulationState"],[3,"String"],[4,"Result"],[15,"str"],[3,"MetaOptions"],[3,"SolarOptions"],[3,"SimulationStateHeader"],[15,"usize"],[3,"TypeId"],[3,"CurrentWeather"],[3,"IRViewFactorSet"],[3,"OpticalInfo"],[3,"Formatter"],[6,"Result"],[3,"SolarSurface"],[3,"Ray3D"],[3,"Vec"],[3,"Scene"],[3,"DCFactory"],[15,"bool"],[6,"Matrix"],[6,"Float"],[8,"Fn"],[3,"Polygon3D"]]},\
"simple_light":{"doc":"","t":[5],"n":["main"],"q":["simple_light"],"d":[""],"i":[0],"f":[[[]]],"p":[]}\
}');
if (typeof window !== 'undefined' && window.initSearch) {window.initSearch(searchIndex)};
if (typeof exports !== 'undefined') {exports.searchIndex = searchIndex};
