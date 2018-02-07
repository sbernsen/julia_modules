


module radarIO

export gps_xmlread, h5co, h5cmp


using HDF5, DataFrames

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
  gps_mat = []

  # Open the hdf5 file
  fid = h5open(h5filename)
  
  if !isempty( fid[line_number ] )
  
    # Allocate space 
    M = length( read(fid[line_number*"/location_0/datacapture_0/echogram_0"]) )
  
    location_id = names(fid[line_number])
    
    N = length(location_id)
    
    # Allocate Space 
    trace_mat = zeros(M, N)
    
    
    for j = 1:N
    
        attr_name = line_number*"/"*location_id[j]*"/datacapture_0/echogram_0"
        trace_h5 = fid[attr_name]
        trace_mat[:,j] = trace_mat[:,i] + read(trace_h5)

        gps_mat = push!(gps_mat, gps_xmlread( read(trace_h5[gps_name]) ) )

    end
    
    gps_mat = DataFrame(gps_mat, names = location_id)
  end

  return trace_mat, gps_mat
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
                trace_mat[:,i] = trace_mat[:,i] + read(trace_h5)
                
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


end

