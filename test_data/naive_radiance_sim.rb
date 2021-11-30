#!/usr/bin/ruby

@INCLUDE_SUN = false
@CORES = 4

@WEAFILE = "./Wellington.wea"
@SENSORFILE = "./sensors.pts"
@N_HOURS = 24 * 7

# Read WEA file
weather_data = File.readlines(@WEAFILE);
@PLACE = weather_data.shift.strip.split(" ")[-1]
@LATITUDE = weather_data.shift.strip.split(" ")[-1].to_f
@LONGITUDE = weather_data.shift.strip.split(" ")[-1].to_f
@STANDARD_MER = weather_data.shift.strip.split(" ")[-1].to_f
@ELEVATION = weather_data.shift.strip.split(" ")[-1].to_f
@UNITS = weather_data.shift.strip.split(" ")[-1].to_i

# Read Sensors
sensors = File.readlines(@SENSORFILE)
@N_SENSORS = sensors.length

weather_data.each_with_index{|ln, index|    
    break if index >= @N_HOURS

    month, day, hour, direct_normal, diffuse_horizontal = ln.strip.split(" ").map{|x| x.to_f}
    month = month.to_i
    day = day.to_i    

    if diffuse_horizontal > 0 then
        # day time

        # Write sky
        sky_filename = "./perez.sky"
        File.open(sky_filename, 'w'){|f|
            f.puts "!gendaylit #{month} #{day} #{hour} #{@INCLUDE_SUN ? "" : "-s" } -g 0 -a #{@LATITUDE} -o #{@LONGITUDE} -m #{@STANDARD_MER} -i 60 -W #{direct_normal} #{diffuse_horizontal}"

            f.puts "skyfunc glow skyglow
            0
            0
            4 1 1 1 0
            
            skyglow source sky
            0
            0
            4 0 0 1 180
            
            skyfunc glow grndglow
            0
            0
            4 1 1 1 0
            
            grndglow source ground
            0
            0
            4 0 0 -1 180"
        }

        # Create octree
        octree_name = "octree.oct"
        `oconv #{sky_filename} room.rad > #{octree_name}`

        # rtrace
        results_file = "temp_results.tmp"
        ambient_file = "ambient.amb"
        `cat #{@SENSORFILE} | rtrace -h -n #{@CORES} -af #{ambient_file} -I+ -ab 8 -ad 6400 #{octree_name} > #{results_file}`
        results = File.readlines(results_file).map{|r|
            r.split(" ")[0].to_f * 179
        }
        puts results.join(",")
        

        # Clean!
        File.delete(sky_filename)
        File.delete(octree_name)
        File.delete(results_file)
        File.delete(ambient_file)
        
    else
        # night time
        puts Array.new(@N_SENSORS, 0.0).join(",")
    end

    
}