module MainDriver

using DrWatson
using TimerOutputs

@quickactivate "GridapODEsTests.jl"

# Here you may include files from the source directory
include(srcdir("beltrami_flow.jl"))

println(
"""
Currently active project is: $(projectname())
Path of active project: $(projectdir())

"""
)

using .BeltramiFlow

# Initial run
params = BeltramiParams(vtk_output=true)
beltrami_flow_solver(params)

# Timings
const to = TimerOutput()

# ODE solvers
ode_solver_types = [:ThetaMethod]#,:GeneralizedAlpha,:RK_CN,:RK_SDIRK]
dts = [0.01/(2^i) for i in 0:3]

eus = Float64[]
eps = Float64[]
for ode_solver_type in ode_solver_types
  for dt in dts
    case = "$ode_solver_type"*"-$dt"
    params = BeltramiParams(dt=dt, ode_solver_type=ode_solver_type, ne=4)
    eᵤ,eₚ = @timeit to case beltrami_flow_solver(params)
    push!(eus,eᵤ)
    push!(eps,eₚ)
  end
end

println("eᵤ = $eus")
println("eₚ = $eps")

end
