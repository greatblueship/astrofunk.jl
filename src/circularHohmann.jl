function circHohmann(r₁,r₂,μ)
#[dV1,dV2,T] = circHohmann(r1,r2,μ)
#   calculates circular Hohmann transfer details given two radii
#   input and output are the same units
#
# Alex Hoffman
# 10/31/2020

aₜ   = 0.5(r₁ + r₂) #transfer semi-major axis
v1 = sqrt(μ/r₁)
v2 = sqrt(μ/r₂)
vt1 = sqrt(μ*(2/r₁-1/aₜ))
vt2 = sqrt(μ*(2/r₂-1/aₜ))
dV1 = vt1 - v1
dV2 = v2 - vt2
T = 2pi*sqrt(aₜ^3/μ)

return [dV1, dV2, T]
end
