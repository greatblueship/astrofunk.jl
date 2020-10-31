ang = pi/3
c = DCM(ang,'x',false)
b = [6;2;6]
d = c*c*b
print(d)

Earth_a = 149600000
Mars_a = 227920000
Sun_mu = 1.3271*10^11
e = circHohmann(Earth_a, Mars_a, Sun_mu)
