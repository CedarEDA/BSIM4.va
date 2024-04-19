"""
This module contains `bsim4_va`, the path to BSIM4 in Verilog-A for use in simulators.
It also exports `bsim4`, which is the Julia lowering of the BSIM4 model.
"""
module BSIM4

using RelocatableFolders 
using CedarSim

# This will hold the path to our file
const bsim4_va = @path joinpath(@__DIR__, "bsim4.va")

Base.include(@__MODULE__, VAFile(bsim4_va))

export bsim4

end # module
