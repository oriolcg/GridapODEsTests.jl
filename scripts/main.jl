module MainDriver

using DrWatson
using TimerOutputs
using Plots

@quickactivate "GridapODEsTests.jl"

# Here you may include files from the source directory
include(srcdir("beltrami_flow.jl"))
include(srcdir("manufactured_linear.jl"))

println(
"""
Currently active project is: $(projectname())
Path of active project: $(projectdir())

"""
)

using .BeltramiFlow
using .ManufacturedLinear

# Initial run
params = ManufacturedLinearParams(vtk_output=true)
manufactured_linear_solver(params)

# Timings
const to = TimerOutput()

# ODE solvers
ode_solver_types = [:ThetaMethod]#,:GeneralizedAlpha,:RK_CN,:RK_SDIRK]
dts = [0.1/(2^i) for i in 0:3]

for ode_solver_type in ode_solver_types

  # Initialize output vectors
  eus = Float64[]
  eps = Float64[]
  times = Float64[]
  memory = Float64[]

  # Time convergence
  for dt in dts
    case = "$ode_solver_type"*"-$dt"
    params = ManufacturedLinearParams(dt=dt, ode_solver_type=ode_solver_type, ne=4)
    eᵤ,eₚ = @timeit to case manufactured_linear_solver(params)
    push!(eus,eᵤ)
    push!(eps,eₚ)
    println(to)
    push!(times,TimerOutputs.time(to[case]))
    push!(memory,TimerOutputs.allocated(to[case]))
  end

  eu_vs_dt = plot(dts,eus,xaxis=:log10,yaxis=:log10,label="$ode_solver_type")
  plot!(eu_vs_dt,dts,1.0e-2 .* dts .^2,label="dt²",color=:black,ls=:dash)
  plot!(eu_vs_dt,dts,1.0e-1 .* dts .^3,label="dt³",color=:black,ls=:dashdot)
  eu_vs_dt_file = plotsdir("eu_vs_dt-$ode_solver_type")
  savefig(eu_vs_dt,eu_vs_dt_file)
end

end
