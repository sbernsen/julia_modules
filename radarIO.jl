


module radarIO

export read_steps, preprocessing_steps, get_lines, gps_xmlread, h5co, h5cmp, bsdeg2decdeg



using HDF5, DataFrames, Preprocessing

##### -------------------------------------------------------------------- #####
#                               Object Definitions                             #
##### -------------------------------------------------------------------- #####

# Define the object returned in read_steps()
struct preprocessing_steps
    filter_type::String
    fc::Any
    fs::Float64
    d_rxtx::Float64
    dxi::Float64
    rms_norm::Bool
    stat_norm::Bool
    kwindow::Int64
    stat_type::String
    rm_bgrnd::Bool
    stack::Bool
    stack_window::Int64
end



function read_steps(survey_type)
# Read the steps for preprocessing. The function looks for a file names 'preprocessing_steps.txt'
#
# survey_type - either common offset "co", common mid/depth-point "co"
#
# -----------------------------------------------------------------------------#

    pre_steps = readdlm("preprocessing_steps.txt", ':')

    rms_norm = (pre_steps[ pre_steps[:,1] .== "rms_norm", 2])[1]
    stat_norm = (pre_steps[ pre_steps[:,1] .== "stat_norm", 2])[1] 

    if stat_norm 
        stat_type = strip( (pre_steps[ pre_steps[:,1] .== "stat_type", 2] )[1] )
        kwindow = (pre_steps[ pre_steps[:, 1].== "kwindow", 2])[1]
    else 
        kwindow = 0 
        stat_type = "none"
    end

    rm_bgrnd = (pre_steps[ pre_steps[:,1] .== "rm_bgrnd", 2])[1]
    
    # We need the cmp or wide angle reciever interval or the common offset 
    d_rxtx = (pre_steps[ pre_steps[:,1] .== "d_rxtx", 2] )[1]

    fs = (pre_steps[ pre_steps[:,1] .== "fs", 2])[1] 
    filter_type = strip( (pre_steps[ pre_steps[:,1] .== "filter_type", 2])[1] )
    

    if filter_type == "bp"
        fc = strip( ( pre_steps[ pre_steps[:,1].== "corner", 2] )[1] )
        fc = [ float( split(fc, '/' )[1] ), float( split(fc, '/' )[2]) ]./fs
    elseif filter_type == "none"
        fc = 1.0
    else
        fc = (pre_steps[ pre_steps[:,1] .== "corner", 2])[1]
        fc = float(fc)/fs
    end

    if survey_type != "co"
      # get the initial offset of the reciever and transmitter
      dxi = (pre_steps[ pre_steps[:,1] .== "dxi", 2] )[1]
      stack = false
      stack_window = 1
      return_obj = preprocessing_steps(filter_type, fc, fs, d_rxtx, dxi, rms_norm, stat_norm, kwindow, stat_type, rm_bgrnd, stack, stack_window)
    else   
      stack = (pre_steps[ pre_steps[:,1] .== "stack", 2])[1]
      stack_window = (pre_steps[ pre_steps[:,1] .== "stack_window", 2])[1]
      return_obj = preprocessing_steps(filter_type, fc, fs, d_rxtx, 0.0, rms_norm, stat_norm, kwindow, stat_type, rm_bgrnd, stack, stack_window)
    end

    return return_obj
end



function get_lines(filename)
  fid = h5open(filename) 
  lines = names(fid)

  return lines 

end


function gps_xmlread(metadata_xml) 

    no_sat = "Num _Sat" 
    time_stamp = "GPS_timestamp_UTC" 
    lat = "Lat_N"
    lon = "Long_ W" 
    quality = "Fix_Quality" 
    ele = "Alt_asl_m"
    gps_fix = "GPS Fix valid"
    gps_msg = "GPS Message ok"
    
    
    metadata = split(metadata_xml, r"<Cluster\r\n|</Cluster>\r\n|<Name>|</Name>|\r\n|<String>|</String>|<Boolean>|</Boolean>", keep = false )
    
    
    no_sat = split( metadata[ find(metadata .== no_sat)+1][1], r"<Val>|</Val>", keep = false)
    
    if !isempty(no_sat)
        no_sat = no_sat[1]
    else
        no_sat = "NA"
    end
    
    # GPS time stamp
    time_stamp = split( metadata[ find(metadata .== time_stamp)+1][1], r"<Val>|</Val>", keep = false)
    
    if !isempty(time_stamp)
        time_stamp = time_stamp[1]
    else
        time_stamp = "NA"
    end
    
    # A quality indicator of the gps
    quality = split( metadata[ find(metadata .== quality)+1][1], r"<Val>|</Val>", keep = false)
    
    if !isempty(quality)
        quality = quality[1]
    else
        quality = "NA"
    end
    
    # whether or not there was a satellite fix
    gps_fix = split( metadata[ find(metadata .== gps_fix)+1][1], r"<Val>|</Val>", keep = false)
    
    if !isempty(gps_fix)
        gps_fix = gps_fix[1]
    else
        gps_fix = "NA"
    end
    
    # GPS messages that indicate usable data
    gps_msg = split( metadata[ find(metadata .== gps_msg)+1][1], r"<Val>|</Val>", keep = false)
    
   if !isempty(gps_msg)
        gps_msg = gps_msg[1]
    else
        gps_msg = "NA"
    end 
    
    !!!!
    
    # Latitude
    lat = split( metadata[ find(metadata .== lat)+1][1], r"<Val>|</Val>", keep = false)
    
    if !isempty(lat)
        lat = lat[1]
    else
        lat = "NA"
    end
    
    
    # Longitude 
    lon = split( metadata[ find(metadata .== lon)+1][1], r"<Val>|</Val>", keep = false)
    if !isempty(lon)
        lon = lon[1]
    else
        lon = "NA"
    end
    
    # Elevation 
    ele = split( metadata[ find(metadata .== ele)+1][1], r"<Val>|</Val>", keep = false)
    
    if !isempty(ele)
        ele = ele[1]
    else
        ele = "NA"
    end
    
    metadata = [no_sat, time_stamp, quality, gps_fix, gps_msg, lat, lon, ele]
    
    return metadata
end


function h5co(h5filename, line_number)

  # Input line_number as "line_##" where ## is the number (e.g. "line_1" or "line_21")
  
  # Define constants
  file_root = split(h5filename, ".h5")[1]
  gps_root = file_root*"_gps"
  gps_name = "GPS Cluster- MetaData_xml"
  gps_header = ["no_sat", "time_stamp", "quality", "gps_fix", "gps_msg", "lat", "lon", "ele"]

  # Open the hdf5 file
  fid = h5open(h5filename)
  
  if !isempty( fid[line_number ] )
  
    # Allocate space 
    M = length( read(fid[line_number*"/location_0/datacapture_0/echogram_0"]) )
  
    location_id = names(fid[line_number])
    
    N = length(location_id)
    
    # Allocate Space 
    trace_mat = zeros(M, N)
    gps_mat = Array{String}(N, length(gps_header) )
    
    for j = 1:N
    
        attr_name = line_number*"/"*location_id[j]*"/datacapture_0/echogram_0"
        trace_h5 = fid[attr_name]
        trace_mat[:,j] = read(trace_h5)

        gps_mat[j,:] = gps_xmlread( read(trace_h5[gps_name]) )

    end
    
    ind = find(  gps_mat[:,4] .!= "0")
    
    if !isempty(ind)
      trace_mat = trace_mat[:,ind] 
      gps_mat = gps_mat[ind,:]
      location_id = location_id[ind] 
    end
    
    ind = find( gps_mat[:,5] .!= "0")
    if !isempty(ind)
      trace_mat = trace_mat[:,ind] 
      gps_mat = gps_mat[ind,:]
      location_id = location_id[ind] 
    end

  end

  return trace_mat, gps_mat, location_id
end


function h5cmp(h5filename)
# Each line consists of a set of time series for an antenna location that needs to be stacked. This is 

  gps_header = ["no_sat", "time_stamp", "quality", "gps_fix", "gps_msg", "lat", "lon", "ele"]
  gps_name = "GPS Cluster- MetaData_xml"

  fid = h5open(h5filename)
  line_id = names(fid)
  
  # Allocate space 
  M = length( read(fid[line_id[2]*"/location_0/datacapture_0/echogram_0"]) )
  n = length(line_id) 
  
  trace_mat = zeros(M, n)
  lle = zeros(n, 3 )
  keep_cols = []
  
  for i in 1:n

    if !isempty( fid[line_id[i] ] ) 
        keep_cols = push!(keep_cols, i) 

        location_id = names(fid[line_id[i]])
        m = length(location_id)
        gps_mat = []
        
        for j in 1:m
                attr_name = line_id[i]*"/"*location_id[j]*"/datacapture_0/echogram_0"
                trace_h5 = fid[attr_name]
                temp_trace = read(trace_h5)

                if Int(length(temp_trace)/size(trace_mat, 1) ) > 1
                    temp_trace = decimate(temp_trace, Int(length(temp_trace)/size(trace_mat, 1) ) )
                end

                trace_mat[:,i] = trace_mat[:,i] + temp_trace
                
                gps_mat = push!(gps_mat, gps_xmlread( read(trace_h5[gps_name]) ) )
        end
        
        gps_mat = DataFrame(gps_mat)
         
        #Check rows 4 and 5 to see if data points are valid
        gps_mat = gps_mat[:, find( convert(Array, gps_mat[4,:]) .!= "0") ]
        gps_mat = gps_mat[:, find( convert(Array, gps_mat[5,:]) .!= "0") ]

        lle[i,:] = mapslices( mean, float(convert( Array, gps_mat[6:8,:] ) ), 2) 

                
        trace_mat[:,i] = trace_mat[:,i]./m
        
    end
  
  end
  
  trace_mat = trace_mat[:,keep_cols] 
  lle = lle[keep_cols,:]

  return trace_mat, lle
end



function bsdeg2decdeg(lat, lon)
#= 
The blue ice systems radar provides gps coordinates in the format of degmm.mmmm
for latitude and longitude. We need to convert it to a usable format. 

!!!!! Longitude is always going to be returned as postive regardless of East/West !!!!!
=#

  deg_lat = round.(lat, -2) 
  deg_lon = round.(lon, -2) 

  min_lat = (lat - deg_lat)./60
  min_lon = (lon - deg_lon)./60

  deg_lat = deg_lat./100 
  deg_lon = deg_lon./100 

  lat = deg_lat + min_lat 
  lon = deg_lon + min_lon 

  return lat, lon


end



# Like all things, the module must come to an end
end