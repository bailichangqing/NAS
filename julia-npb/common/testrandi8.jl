
include("./randi8.jl")
import Randi8


x = Randi8.randlc(1.0,1.0)


for i in 1:100
	x = Randi8.randlc(x, x)
	@printf("Randi8[%d] = %f\n", i, x)
end
