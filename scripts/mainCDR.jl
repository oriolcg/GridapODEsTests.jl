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
uₕ = periodic_convection_diffusion_solver(params)
params = PeriodicConvectionDiffusionParams(ne=3,Δt=0.001)
uₕ_ref = periodic_convection_diffusion_solver(params)

e = abs(uₕ-uₕ_ref)
println("Finished successfully with error: $e")

end
