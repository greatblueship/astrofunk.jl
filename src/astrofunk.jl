# Alex Hoffman
module astrofunk

# include the file and then export the functions
include("circularHohmann.jl")
export circHohmann

# include("dirCosMat.jl")
# export DCM
#
# include("jacobiConstant.jl")
# export JC
# export xtest

include("basicEquations.jl")
export JC
export DCM
export CR3BP_EOM
export U_star

end #astrofunk module
