# module providing various methods for Frequency Time Analysis

module FTAN 

# --------------------------------------------------------------------------- #

export cwt, invcwt, wave_bases, stran_hyp, stran_gaus

# --------------------------------------------------------------------------- #

using ToeplitzMatrices, ImageFiltering, RCall

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Code was adopted from matlab code written by Christopher Torrence and Gilbert
# P. Compo. Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
# Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.

function wave_bases(mother, k, sc, param)
#=
Generate the daughter wavelet given the scale, wavenumber, and type specific
parameter. See Table 1 in the above ref. for equations of the Morlet and DoG
wavelets

Input Variables:
	mother - string value for the mother wavelet 
	k - the wavenumber 
	sc - the scale 
	param - k0 for morlet; m (derivative) for dog 

Output Variables:
	daughter - the daughter wavelet 

=#

	n = length(k)

	if mother == "morlet"

		expnt = -( (sc.*k - param).^2)/2.*Cint.(k.>0)
		nrm = sqrt( sc*k[2])*(pi^-0.25)*sqrt(n)
		daughter = expnt.*exp.(expnt).*Cint.(k.>0)

	elseif mother == "dog"

		expnt = -(sc.*k).*Cint.(k.>0)
		nrm = sqrt(sc*k[2])*( (2^param)/sqrt(param*factorial(2param-1) ) )*sqrt(n) 
		daughter = (nrm*( (sc.*k).^param) ).*Cint.(k.>0)

	end 

	return daughter

end 


# --------------------------------------------------------------------------- #

function cwt(x, dsc, nsc, mother, param)
#=
Compute the continous wavelet transform for a normalized period. Multiply 
period by dt to get the true period of the time series
Input Variables:
	x - 
	dt - 
	dsc - scale spacing; smaller d_sc gives better resolution but more time to plot
	nsc - the number of scales 
	mother - the name of the mother wavelet

=#
	n = length(x) 
	dt = 1

	# Create the scale array
	if isempty(dsc)
		dsc = 1/4 # scale spacing
	end

	
	NSC =  Int(ceil( (log2(n/2) )/dsc ) )

	if isempty(nsc)
		nsc = NSC 
	elseif nsc > NSC 
		nsc = NSC 
	end

	# remove the mean from the time series
	x = x - mean(x) 

	# pad the time series 
	x = [zeros(n, 1); x]

	# compute fft 
	X = fft(x) 

	# create the scale array
	sc = collect( 2*(2*dt).^( (0:nsc)*dsc) ) 
	
	# create the wavenumber array
	k = collect(1:n).*( (pi/n*dt) )
	k = [0; k; -k[ Int(floor( (2n-1)/2 ) ):-1:1] ]

	# define the fourier factor, cone of influence, and degrees of freedom.
	if mother == "morlet"
		if isempty(param)
			param = 6 
		end
		fourier_factor = (4*pi)/(param + sqrt(2 + param^2) )
		coi = (fourier_factor/sqrt(2) )*dt*[1e-5; collect(1:( (n+1)/2 - 1) ); collect( ( (n+1)/2 - 1):-1:1); 1e-5]
		dofmin = 2	
	elseif mother == "dog"

		if isempty(param) 
			param = 2
		end

		fourier_factor = 2pi*sqrt(2/(2*param+1) )
		coi = (fourier_factor/sqrt(2) )*dt*[1e-5; collect(1:( (n+1)/2 - 1) ); collect( ( (n+1)/2 - 1):-1:1); 1e-5]
		dofmin = 1
	else
		error("Mother wavelet must be specified as 'morlet' or 'dog' ")
	end


	# allocate space 
	period = zeros( size(sc) )
	wave = zeros(nsc+1, 2n)

	for i = 1:(nsc+1)
		daughter = wave_bases(mother,k,sc[i],param)
		wave[i,:] = real(ifft(X.*daughter))
	end

	wave = wave[:,(n+1):end]
	period = fourier_factor.*sc;

	return wave, period, sc, coi, nsc, k
end


# --------------------------------------------------------------------------- #

function invcwt(wave, mother, sc, param, k)
# Compute the inverse cwt. There's some bad coding that I rewrote so I need to check it and modify
	# Get the real part of the wavelet transform 
	Wr = real(wave) 
	n = size(Wr, 2) 

	s = repmat(sc, [1, size(wave, 2) ])
	summand = sum( Wr./sqrt(s), 1)
	Wdelta = zeros( size(sc) )

	for i = 1:length(sc)
		daughter = wave_bases(mother, k, sc[i], param) 
		Wdelta[i] = (1/N)*sum(daughter) 
	end 

	C = sum( real(Wdelta)./sqrt(scale) )
	Xrec = (1/C)*summand

	return Xrec

end


# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #


function stran_hyp(ts, gamma_f, gamma_b, lambda)
# Compute S-transform using a hyperbolic gaussian window. See Pinnegar and Mansinha 2003 for details.
#
# Input Variables:
#	ts
#	gamma_f - forward taper paramater
#	gamma_b - backward taper parameter
#	lambda - positive curvature parameter 
#
# Output Variables:
#
#

	M = length(ts)
	H = [fft(ts); fft(ts) ]
	
	# Define positive and negative time (t - τ)
	t = ([ 1:floor(M/2); -floor(M/2)+1:0]) - (gamma_b-gamma_f)/2*sqrt(lambda/gamma_f/gamma_b)
	
	X = (gamma_b+gamma_f)/(2*gamma_b*gamma_f)*t + (gamma_b-gamma_f)/(2*gamma_b*gamma_f)*sqrt.(t.^2 + lambda) 
	
	W = [1; zeros(M-1,1)]


	f = collect(0:(floor(M/2)-1) )
	s_tran = zeros(length(f),M)

	for i = 2:length(f)
		s_tran[i,:] = real(ifft(H[i:(i+M-1)].*W))
		w = (2*abs(f[i])/(sqrt(2pi)*(gamma_b+gamma_f) ) )*exp.(-( (f[i]^2)*X.^2)/(2*M^2) )
		W = fft( w/sum(w) )
	end

	return f, s_tran
end

# --------------------------------------------------------------------------- #

function stran_gaus(ts)
# Compute S-Transform  to a normalized frequency

	# ts is a 1xN one-dimensional series
	N = length(ts)
	n=Int(floor(N/2))

	odvn = Int(N-2n) # zero if N is positive and 1 if N is negative 

	f=[ collect(0:n); collect(-n+1-odvn:-1)]'./N

	H=fft(ts)
	invfk=(1./f[2:n+1] )'

	#Compute all frequency domain Gaussians as one matrix
	W=2*pi*repmat(f,n,1).*repmat(invfk,1,N)

	#H Gaussian in freq domain
	G=exp.((-W.^2)/2); 


	## End of frequency domain Gaussian computation

	#Compute Toeplitz matrix with the shifted fft(ts)
	HW=Toeplitz(H[1:n+1],H)

	#Exclude the first row, corresponding to zero frequency
	HW=HW[2:n+1,:]

	#Compute Stockwell Transform
	s_tran=real( ifft(HW.*G,2) )

	#Add the zero freq row
	st0=mean(ts)*ones(1,N)

	s_tran=[st0;s_tran]
	f = [0; 1./invfk]

	return f, s_tran
end



# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #




end


