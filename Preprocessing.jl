module Preprocessing


export read_steps, preprocessing_steps, sd, remove_redundancy, running_stat, runsd_normalize, background_removal, rms_norm, butt_design

##### -------------------------------------------------------------------- #####
#                            PREPROCESSING_FUNCTIONS                           #
#
# Created by Steven Bernsen with hopes of someday graduating
# University of Maine
#
# This contains simple preprocessing functions to filter, remove noise, 
# normalize, etc.
##### ------------------------------------------------------------------- #####

using Indicators, RCall, DSP, CSV


##### -------------------------------------------------------------------- #####
#                               Object Definitions                             #
##### -------------------------------------------------------------------- #####

# Define the object returned in read_steps()
struct preprocessing_steps
    filter_type::String
    fc::Array
    fs::Float64
    d_rxtx::Float64
    rms_norm::Bool
    stat_norm::Bool
    kwindow::Int64
    stat_type::String
    rm_bgrnd::Bool
end



function read_steps()
# Read the steps for preprocessing. The function looks for a file names 'preprocessing_steps.txt'
#
# 
# -----------------------------------------------------------------------------#

pre_steps = CSV.read("preprocessing_steps.txt", header = [], delim = ':')

rms_norm = parse( Bool, strip( (pre_steps[ pre_steps[:,1] .== "rms_norm", 2])[1] ) )
stat_norm = parse( Bool, strip( (pre_steps[ pre_steps[:,1] .== "stat_norm", 2])[1] ) )

if stat_norm 
    stat_type = strip( (pre_steps[ pre_steps[:,1] .== "stat_type", 2] )[1] )
    kwindow = parse( Int, strip( (pre_steps[ pre_steps[:, 1].== "kwindow", 2])[1] ) )
else 
    kwindow = "NA" 
    stat_type = "NA"
end

rm_bgrnd = parse( Bool, strip( (pre_steps[ pre_steps[:,1] .== "rm_bgrnd", 2])[1] ) )
d_rxtx = float( strip( (pre_steps[ pre_steps[:,1] .== "d_rxtx", 2] )[1] ) )
fs = float( strip( (pre_steps[ pre_steps[:,1] .== "fs", 2])[1] ) )
filter_type = strip( (pre_steps[ pre_steps[:,1] .== "filter_type", 2])[1] )
fc = strip( ( pre_steps[ pre_steps[:,1].== "corner", 2] )[1] )

if filter_type == "bp"

    fc = [ float( split(fc, '/' )[1] ), float( split(fc, '/' )[2]) ]./fs
elseif filter_type == "none"
    fc = 1.0
else
    fc = float(fc[1])/fs
end

return_obj = preprocessing_steps(filter_type, fc, fs, d_rxtx, rms_norm, stat_norm, kwindow, stat_type, rm_bgrnd)
return return_obj

end



function sd(v)
# sd computes the standard deviation of a set of numbers. The intrinsic function
# std won't work. Most of the time I really like julia but not when trying to
# compute the standard deviation or variance of an array of type any
#
# --------------------------------------------------------------------------- #

sd = sqrt( sum( (v - mean(v) ).^2 )/(k-1) )
return sd

end



function remove_redundancy(data_array, threshold)
# remove_redundancy removes the locations that are repeats using correlation
#
# Input Variables:
#   data_array - the matrix/data frame of values to process
#
# Output Variables:
#
# ---------------------------------------------------------------------------- #

n = size(data_array, 2)
corr_value = zeros(n-1, 1)

for i in 2:n
    corr_value[i-1] = cor(data_array[:,i-1], data_array[:,i])
end

data_clean = data_array[:, find(corr_value .<= threshold) ]
return data_clean, corr_value

end



function running_stat(ts, k, stat_type)
# running_stat computes the running statistic of a series for a specified window
#
# Input Variables:
#   ts - the m-by-1 timeseries to be normalized
#   k - the window length
#   stat_type - the type of statistic to be computed; "mean", "std"
#https://www.google.com/search?client=ubuntu&channel=fs&q=gray+scale+image+plot+R&ie=utf-8&oe=utf-8
# Output Variables:
#   run_stat - the running statistic
#
# ---------------------------------------------------------------------------- #

# Allocate space
n = length(ts)
run_stat = zeros(n,1 )

# the window length must be odd
if iseven(k)
    k2 = Int( k/2 )
    k = Int( k + 1 )
else
    k2 = Int( floor(k/2) )
end

# Compute the running mean and standard deviation at the ends
if stat_type == "mean"
    for i = 1:k2
        run_stat[i] = mean( ts[1:(i+k2)] )
    end
    for i = (n-k2+1):n
        run_stat[i] = mean( ts[(i-k2):n] )
    end
    # compute for the full window length
    for i = (1+k2):(n-k2)
        run_stat[i] = mean( ts[(i-k2):(i+k2)] )
    end
elseif stat_type == "abs_mean"
    ts_abs = abs(ts)
    for i = 1:k2
        run_stat[i] = mean( ts_abs[1:(i+k2)] )
    end
    for i = (n-k2+1):n
        run_stat[i] = mean( ts_abs[(i-k2):n] )
    end
    # compute for the full window length
    for i = (1+k2):(n-k2)
        run_stat[i] = mean( ts_abs[(i-k2):(i+k2)] )
    end
else
    for i = 1:k2
        run_stat[i] = std( ts[1:(i+k2)] )
    end
    for i = (n-k2+1):n
        run_stat[i] = std( ts[(i-k2):n] )
    end
    # compute for the full window length
    for i = (1+k2):(n-k2)
        run_stat[i] = std( ts[(i-k2):(i+k2)] )
    end
end

ts_norm = ts./run_stat

return ts_norm


end



function background_removal(data_array)
# background_removal computes a mean stack to estimate the background noise or
# ringing in the timeseries. The noise is then removed from the array of data
# and returned.
#
# Input Variables:
#
# Output Variables:
#
#
# ---------------------------------------------------------------------------- #

## lets get the mean values
n = size(data_array,2)

# first let's scale each timeseries
#for i = 1:n
#    data_array[:,i] = ( data_array[:,i] - mean(data_array[:,i]) )./std(data_array[:,i] )
#end

# sum all of the timeseries, there isn't a function that supports data frames
data_avg = data_array[:,1]
for i = 2:n
    data_avg = data_avg + data_array[:,i]
end
data_avg = data_avg./n

# difference the original data with the average value
for i = 1:n
    data_array[:,i] = data_array[:,i] - data_avg
end

return data_array, data_avg

end



function rms_norm(ts)
# Compute the cross correlation in the frequency domain of each column
#
# Input Variables:
#
# Output Variables:
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n = length(ts)
rms = sqrt( (1/n)*(ts'*ts) )

ts = ts./rms

return ts

end



##### -------------------------------------------------------------------- #####
#                              FILTERING_FUNCTIONS                             #
#
# Created by Steven Bernsen because the wheel still isn't good enough
# University of Maine
##### -------------------------------------------------------------------- #####


function butt_design( filter_type, fc )
# lowpass_butter applies a low pass, 2 pole butterworth filter to a time series.
# The filter is designed to a normalized sampling frequency of 1 and can be
# scaled afterword. For example, the cutoff frequency of 20 MHz for a
# timeseries sampled at 250 MHz would be fc = 20/250 = 0.08
#
# Input Variables:
#   ts - the input timeseries
#   fc - the corner frequency or frequencies. For bandpass filtering, input as
#       [hp_corner, lp_corner]
#   type - specified as "hp", "lp", "bp" for highpass, lowpass and bandpass,
#       respectively
#
# Output Variables:
#
#
##### -------------------------------------------------------------------- #####

if filter_type == "bp"
    responsetype = Bandpass( fc[1], fc[2] )
elseif filter_type == "lp"
    responsetype = Lowpass(fc)
else
    responsetype = Highpass(fc)
end

designmethod = Butterworth(2)
filter_object = digitalfilter(responsetype, designmethod)

return filter_object

end



### Thats it, son!
end



\