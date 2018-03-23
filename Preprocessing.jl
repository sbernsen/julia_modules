
module Preprocessing


export sd, remove_redundancy, running_stat, runsd_normalize, background_removal, 
        rms_norm, butt_design, time_adjust, decimate, whiten, trace_stack

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





function sd(v)
# sd computes the standard deviation of a set of numbers. The intrinsic function
# std won't work. Most of the time I really like julia but not when trying to
# compute the standard deviation or variance of an array of type any
#
# --------------------------------------------------------------------------- #

    sd = sqrt( sum( (v - mean(v) ).^2 )/(k-1) )
    re
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
      ts_abs = abs.(ts)
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


  function time_adjust(ts, dx, fs, v)

    t_adj = Int(round(fs*dx/v) )
    add_zeros = zeros( t_adj, 1) 
    m = length(ts) 

    # make sure that ts is m-by-1 
    ts = reshape(ts, m, 1) 

    ts = [add_zeros; ts] 
    ts = ts[1:m] 

    return ts, t_adj

end


function trace_stack(trace_mat, stack_window)


  n = size(trace_mat, 2) 

  if !isodd(stack_window) 
    stack_window = stack_window + 1
  end 

  ind_start = Int( (stack_window+1)/2)

  a = 1
  b = stack_window

  for i = ind_start:(n-ind_start+1)
      trace_mat[:,i] = mapslices(mean, trace_mat[:,a:b], 2)
      a += 1
      b += 1
  end

  return trace_mat
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

###############################################################
## For the following code don't forget to acknowledge Jacob Svensson @ https://gist.github.com/jcbsv . Some code was modified
#################################################

function cheby1(n, r, wp)
# Chebyshev Type I digital filter design.
#
#    b, a = cheby1(n, r, wp)
#
# Designs an nth order lowpass digital Chebyshev filter with
# R decibels of peak-to-peak ripple in the passband.
#
# The function returns the filter coefficients in length
# n+1 vectors b (numerator) and a (denominator).
#
# The passband-edge frequency wp must be 0.0 < wp < 1.0, with
# 1.0 corresponding to half the sample rate.
#
#  Use r=0.5 as a starting point, if you are unsure about choosing r.
  h = digitalfilter(Lowpass(wp), Chebyshev1(n, r))
  tf = convert(PolynomialRatio, h)
  coefb(tf), coefa(tf)
end


function decimate(x, r)
# Decimation reduces the original sampling rate of a sequence
# to a lower rate. It is the opposite of interpolation.
#
# The decimate function lowpass filters the input to guard
# against aliasing and downsamples the result.
#
#   y = decimate(x,r)
#
# Reduces the sampling rate of x, the input signal, by a factor
# of r. The decimated vector, y, is shortened by a factor of r
# so that length(y) = ceil(length(x)/r). By default, decimate
# uses a lowpass Chebyshev Type I IIR filter of order 8.
#
# Sometimes, the specified filter order produces passband
# distortion due to roundoff errors accumulated from the
# convolutions needed to create the transfer function. The filter
# order is automatically reduced when distortion causes the
# magnitude response at the cutoff frequency to differ from the
# ripple by more than 1E–6.

  nfilt = 8
  cutoff = .8 / r
  rip = 0.05  # dB

  function filtmag_db(b, a, f)
    # Find filter's magnitude response in decibels at given frequency.
    nb = length(b)
    na = length(a)
    top = (exp.(-1im*collect(0:nb-1)*pi*f) )'*b
    bot = (exp.(-1im*collect(0:na-1)*pi*f) )'*a
    20*log10(abs(top/bot))
  end

  b, a = cheby1(nfilt, rip, cutoff)
  while all( b == 0) || ( abs( filtmag_db( b, a, cutoff )+rip )>1e-6 )
    nfilt = nfilt - 1
    if nfilt == 0
        break
    end
    b, a = cheby1(nfilt, rip, cutoff)
  end
  y = filtfilt(PolynomialRatio(b, a), x)
  nd = length(x)
  nout = ceil(nd/r)
  nbeg = Int(r - (r * nout - nd))
  y[nbeg:r:nd]
end


##### -------------------------------------------------------------------- #####


function boxcar3(A::AbstractArray)
# Run a 3-by-3-by-... boxcar filter on an array. Courtesy of Tim Holy https://julialang.org/blog/2016/02/iteration
    out = similar(A)
    R = CartesianRange(size(A))
    I1, Iend = first(R), last(R)
    for I in R
        n, s = 0, zero(eltype(out))
        for J in CartesianRange(max(I1, I-I1), min(Iend, I+I1))
            s += A[J]
            n += 1
        end
        out[I] = s/n
    end
    out
end


##### -------------------------------------------------------------------- #####

function whiten(x, freq) 
#=
Spectral whitening of a time series for a given frequency range 

Input Variables:
  x - black Michael Jackson 
  freq - the corner frequencies normalized to 1 unit time sampling frequency to 
      be whitened 
  
Output Variables:
  x - white Michael Jackson
=#
  
  # do a little error checking 
  if sum(Int.(freq .> 0.5) ) != 0  
    error("Seriously? The upper corner is beyond the Nyquist. You need to retake time series.")
  end

  if isempty(freq) 
    freq = [0 0.5]
  end 


  n = length(x) 

  # create the hanning window
  W = hanning(n)
  X = fft(W.*x) 

  mag = abs.(X) 
  mag_mlt = maximum(mag) 

  phase = unwrap(angle.(X))
  f = (collect( 0:(n-1) ) )./n 

  df = abs.(f-freq[1])
  ind1 = find( df .== minimum(df)  )[1]

  df = abs.(f-freq[2])
  ind2 = find( df .== minimum(df)  )[1]

  # Normalize the magnitude across the specified band 
  mag[ ind1:ind2 ] = mag_mlt

  x = real( ifft( mag.*exp.(1im.*phase) ) )

  return x 

end






### Thats it, son!
end