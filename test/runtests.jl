# Alex Hoffman
#cd("/Users/hoffm/.julia/dev")
using astrofunk
using Test

#variables for testing
Earth_a, Mars_a, Sun_mu = 149600000, 227920000, 1.3271*10^11
pos = [-0.270; -0.420; 0] #nondimensional
vel = [0.300; -1.00; 0] #nondimensional
μ₃ = 0.012150664267189

#test loop
@testset "astrofunk.jl" begin
    # Write your tests here.
    @test DCM(pi,'x',false) ≈ [1 0 0; 0 -1 0; 0 0 -1]
    @test circHohmann(Earth_a, Mars_a, Sun_mu) ≈ [2.944004579970198, 2.648345397475325, 44729392.19214477]
    @test JC([pos pos], [vel vel], μ₃) ≈ permutedims([3.186469735193683; 3.186469735193683])
    @test CR3BP_EOM([pos;vel],μ₃,'x') ≈ -0.135541250504792
    @test U_star(pos,μ₃,"yy") ≈ 10.724399149271521
end
