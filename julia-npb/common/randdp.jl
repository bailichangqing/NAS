module NPBRand
# 
# 
#     FUNCTION RANDLC (X, A)
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
#   This routine should produce the same results on any computer with at least
#   48 mantissa bits in double precision floating point data.  On Cray systems,
#   double precision should be disabled.
# 
#   David H. Bailey     October 26, 1990
# 
#      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#      SAVE KS, R23, R46, T23, T46
#      DATA KS/0/
# 
#   If this is the first call to RANDLC, compute R23 = 2 ^ -23, R46 = 2 ^ -46,
#   T23 = 2 ^ 23, and T46 = 2 ^ 46.  These are computed in loops, rather than
#   by merely using the ** operator, in order to insure that the results are
#   exact on all systems.  This code assumes that 0.5D0 is represented exactly.
#


#=================================================================#
#=============           R  A  N  D  L  C             ============#
#=============                                        ============#
#=============    portable random number generator    ============#
#=================================================================#

function randlc(x::Float64, a::Float64)

	const r23 = 0.5 ^ 23
	const r46 = r23 ^ 2
	const t23 = 2.0 ^ 23
	const t46 = t23 ^ 2

#---------------------------------------------------------------------
#   Break A into two parts such that A = 2^23 * A1 + A2.
#---------------------------------------------------------------------

	t1 = r23 * a
	a1 = round(t1)
	a2 = a - t23 * a1
	
#---------------------------------------------------------------------
#   Break X into two parts such that X = 2^23 * X1 + X2, compute
#   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
#   X = 2^23 * Z + A2 * X2  (mod 2^46).
#---------------------------------------------------------------------

	t1 = r23 * x
	x1 = round(t1)
	x2 = x - t23 * x1
	t1 = a1 * x2 + a2 * x1
	t2 = round(r23 * t1)
	z = t1 - t23 * t2
	t3 = t23 * z + a2 * x2
	t4 = round(r46 * t3)
	x = t3 - t46 * t4
	return (r46 * x, x)
end

function vranlc!(n, x::Float64, a::Float64, y::Array{Float64})

	const r23 = 0.5 ^ 23
	const r46 = r23 ^ 2
	const t23 = 2.0 ^ 23
	const t46 = t23 ^ 2

	# Break A into two parts such that A = 2^23 * A1 + A2
	t1 = r23 * a
	a1 = Int(round(t1))
	a2 = a - t23 * a1

	# Generate N results
	for i = 1:n
		t1 = r23 * x
		x1 = Int(round(t1))
		x2 = x - t23 * x1
		t1 = a1 * x2 + a2 * x1
		t2 = Int(round(r23*t1))
		z = t1 - t23 * t2
		t3 = t23 * z + a2 * x2
		t4 = Int(round(r46 * t3))
		x = t3 - t46 * t4
		y[i] = r46 * x
	end

	return x

end

end # !module Randdp!
