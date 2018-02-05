module EGF


export mw_array, array_xcorr, array_xcorr_td, rms_norm

using DSP

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

    # Pad with zeros if needed
    if pad == true
        mw_array = zeros(2*win_len, n)
    end

    for i in 1:n
        mw_array[1:win_len, i] = ts[ a[i]:b[i] ].*win
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
#   arr_xc - the cross correlation of equal dimensions the input arrays
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Compute the FFT on each column
    arr1 = mapslices(fft, arr1, 1)
    arr2 = mapslices(fft, arr2, 1)

    # Compute the dot product
    fft_arr = arr1.*conj(arr2)

    # Get the inverse FFT
    arr_xc = real( mapslices( ifft, fft_arr, 1) )

    # normalize by rms
    arr_xc = mapslices( rms_norm, arr_xc, 1)

    # Stack
    xcf = mapslices(mean, arr_xc, 2)

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
            xcf = xcf + rms_norm( xcorr(ts2[:,i], ts1[:,i] ) )
        end

        xcf = xcf./n

    else
        # Remove the standard deviation
        ts1 = rms_norm(ts1)
        ts2 = rms_norm(ts2)

        xcf = xcorr(ts2[:], ts1[:])
    end


    return xcf

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
