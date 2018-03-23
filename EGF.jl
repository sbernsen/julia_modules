module EGF

export mw_array, array_xcorr, array_xcorr_td, rms_norm

using DSP, Preprocessing

function mw_array(ts, win_len, STEP, pad)
# Assemble an array for the tukey window where each column is a new step in the
# time series then pad the array with zeros if specified.
#
# Input Variables:
#   ts - the m-by-1 time series
#   win_len - the integer value of the window length
#   STEP - the integer step size for the number of
#
# Output Variables:
#
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    # Set the tukey alpha value; 0 - rectangle window, 1 - Hann window
    alpha = 0.2

    # create the tukey window
    win = tukey(win_len, alpha)

    # Lets get the integer values of the time series for each window
    a = [ 1 ]
    b = [ win_len ]
    m = length(ts)
    i = 1

    while (b[i]+STEP) <= m
        a = push!(a, a[i] + STEP )
        b = push!(b, b[i] + STEP )

        i = i + 1
    end

    n = length(a)
    mw_array = zeros(win_len, n)


    for i in 1:n
        mw_array[1:win_len, i] = ts[ a[i]:b[i] ].*win
    end

    # Pad with zeros if needed
    if pad == true
        mw_array = [ zeros(size(mw_array)); mw_array]
    end


    return mw_array

end


function array_xcorr(arr1, arr2)
# Compute the cross correlation in the frequency domain of each column
#
# Input Variables:
#   arr1, arr2 - padded arrays of equal dimensions to be cross correlated
#
# Output Variables:
#   xcf - the cross correlation of equal dimensions the input arrays
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    M = size(arr1, 1) 

    if Bool(mod(M,2))
        m = Int( (M-1)/2 )
    else 
        m = Int( M/2 )
    end

    # Compute the FFT on each column
    arr1 = mapslices(fft, arr1, 1)
    arr2 = mapslices(fft, arr2, 1)

    # Compute the dot product
    fft_causal = arr1.*conj(arr2)
    fft_acausal = arr2.*conj(arr1) 

    # Get the inverse FFT
    arr_c = real( mapslices( ifft, fft_causal, 1) )
    arr_a = real( mapslices( ifft, fft_acausal, 1) )


    # normalize by rms
    arr_c = mapslices( rms_norm, arr_c, 1)
    arr_a = mapslices( rms_norm, arr_a, 1)

    # Stack
    arr_c = ( mapslices(mean, arr_c, 2) )[1:m]
    arr_a = ( mapslices(mean, arr_a, 2) )[m:-1:1]

    xcf = [arr_a; arr_c]

    return xcf

end


function array_xcorr_td(ts1, ts2)
# Compute the cross correlation in the frequency domain of each column
#
# Input Variables:
#   arr1, arr2 - padded arrays of equal dimensions to be cross correlated
#
# Output Variables:
#   arr_xc - the cross correlation of equal dimensions the input arrays
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    m = size(ts1, 1)
    n = size(ts1, 2)

    if n > 1
        xcf = zeros(2m-1, 1)

        for i in 1:n
            xcf = xcf + crosscor(ts2[:], ts1[:], collect(-m:m); demean = true)
        end

        xcf = xcf./n

    else
        xcf = crosscor(ts2[:], ts1[:], collect(-m:m); demean = true)
    end


    return xcf

end


################################################################################
################################# IO Functions #################################
################################################################################



function matrix2SAC(data_arr, sta_list)
#
#
# This function was made specific for 'virtual_gathers.jl' because each gather
# corresponds to the list of stations used to create them.
end


# end the module
end
