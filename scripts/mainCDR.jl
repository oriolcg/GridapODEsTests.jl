module MainCDR

using DrWatson
using TimerOutputs
using Plots

@quickactivate "GridapODEsTests.jl"

# Here you may include files from the source directory
include(srcdir("periodic_convection_diffusion.jl"))

println(
"""
Currently active project is: $(projectname())
Path of active project: $(projectdir())

"""
)

using .PeriodicConvectionDiffusion

# Initial run
params = PeriodicConvectionDiffusionParams(vtk_output=true,ne=3)
periodic_convection_diffusion_solver(params)

end
