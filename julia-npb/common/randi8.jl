
module NPBRand

#---------------------------------------------------------------------
#
#   This routine returns a uniform pseudorandom double precision number in the
#   range (0, 1) by using the linear congruential generator
#
#   x_{k+1} = a x_k  (mod 2^46)
#
#   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
#   before repeating.  The argument A is the same as 'a' in the above formula,
#   and X is the same as x_0.  A and X must be odd double precision integers
#   in the range (1, 2^46).  The returned value RANDLC is normalized to be
#   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
#   the new seed x_1, so that subsequent calls to RANDLC using the same
#   arguments will generate a continuous sequence.
#

function randlc(x::Float64, a::Float64)

	i246m1::Int64 = 0x00003FFFFFFFFFFF
	const d2m46::Float64 = 0.5^46
	
	Lx = Int(round(x))
	La = Int(round(a))

	@printf("Lx is %d\n", Lx)
	Lx = (Lx * La) & i246m1
	randlc = d2m46 * convert(Float64, Lx)
	#@printf("rand returning %f\n", convert(Float64, Lx))
	if convert(Float64, Lx) == 1.0
		@printf("Got 1.0 from initial %f, lx=%d, la=%d, lx*la=%d\n", x, Lx, La, Lx*La)
		
	end
		
	return randlc, convert(Float64, Lx)

end

function vranlc!(n, x, a, y)

	i246m1::Int64 = 0x00003FFFFFFFFFFF
	d2m46  = 0.5^46
	
	Lx = Int(round(x))
	La = Int(round(a))
	
	for i in 1:n
		Lx   = (Lx * La) & i246m1
		y[i] = d2m46 * Float64(Lx)
	end

	return Float64(Lx)

end

end # !module
