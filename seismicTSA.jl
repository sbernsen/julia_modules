
# Compute some common time series analysis on seismic data
using DSP


function psdpdf(ts, window, overlap, fs, ACF)

	# Set the tukey alpha value; 0 - rectangle window, 1 - Hann window
    alpha = 0.2

    # create the tukey window
    taper_win = tukey(window, alpha)

    # Lets get the integer values of the time series for each window
    a = [ 1 ]
    b = [ window ]
    m = length(ts)
    i = 1

    while (b[i]+overlap) <= m
        a = push!(a, a[i] + overlap )
        b = push!(b, b[i] + overlap )

        i = i + 1
    end

    n = length(a)
    window_array = zeros(window, n)


    for i in 1:n
        window_array[1:window, i] = ts[ a[i]:b[i] ].*taper_win
    end

    window_array = [ zeros(size(window_array)); window_array]

    # compute the fft on each column and autocorrelate 
    window_array = mapslices(fft, window_array, 1)

    # We can take the 
    if ACF
    	window_array = ( window_array.*conj(window_array) )[1:window,:]
    	median_psd = [0; mapslices(median, real(window_array), 2) ]
    	std_psd = [0; mapslices(std, real(window_array), 2) ]
   	else 
    	window_array = window_array[1:window,:] 
    	median_psd = [0; mapslices(median, real(window_array), 2) ].^2
    	std_psd = [0; mapslices(std, real(window_array), 2) ]
    end 

    

    freqs = fs.*collect(0:window)./(2*window) 


    return freqs, median_psd, std_psd

end