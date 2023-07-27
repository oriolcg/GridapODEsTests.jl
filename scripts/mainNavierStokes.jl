module MainNavierStokesDriver

using DrWatson
using TimerOutputs
using Plots

@quickactivate "GridapODEsTests.jl"

# Here you may include files from the source directory
include(srcdir("beltrami_flow.jl"))
include(srcdir("manufactured_NavierStokes_linear.jl"))

println(
"""
Currently active project is: $(projectname())
Path of active project: $(projectdir())

"""
)

using .BeltramiFlow
using .ManufacturedNavierStokesLinear

# Initial run
params = ManufacturedNavierStokesLinearParams(vtk_output=true)
manufactured_NavierStokes_linear_solver(params)
params = ManufacturedNavierStokesLinearParams(ode_solver_type=:GeneralizedAlpha)
manufactured_NavierStokes_linear_solver(params)
# params = ManufacturedNavierStokesLinearParams(ode_solver_type=:RK_CN)
# manufactured_NavierStokes_linear_solver(params)
# params = ManufacturedNavierStokesLinearParams(ode_solver_type=:RK_SDIRK)
# manufactured_NavierStokes_linear_solver(params)

# Timings
const to = TimerOutput()

# ODE solvers
ode_solver_types = [:ThetaMethod,:GeneralizedAlpha]#,:RK_CN,:RK_SDIRK]
dts = [0.1/(2^i) for i in 0:3]

# Initialize plots
eu_vs_dt = plot(title="eᵤ vs dt")
eu_vs_time = plot(title="eu vs time")
time_vs_dt = plot(title="time vs dt")
eu_vs_mem = plot(title="eᵤ vs memory")

# Loop over solvers
for ode_solver_type in ode_solver_types

  # Initialize output vectors
  eus = Float64[]
  eps = Float64[]
  times = Float64[]
  memory = Float64[]

  # Time convergence
  for dt in dts
    case = "$ode_solver_type"*"-$dt"
    params = ManufacturedNavierStokesLinearParams(tf=0.2,dt=dt, ode_solver_type=ode_solver_type, ne=20)
    eᵤ,eₚ = @timeit to case manufactured_NavierStokes_linear_solver(params)
    push!(eus,eᵤ)
    push!(eps,eₚ)
    println(to)
    push!(times,TimerOutputs.time(to[case]))
    push!(memory,TimerOutputs.allocated(to[case]))
  end

  # Plotting
  plot!(eu_vs_dt,dts,eus,xaxis=:log10,yaxis=:log10,label="$ode_solver_type")
  plot!(eu_vs_time,times,eus,xaxis=:log10,yaxis=:log10,label="$ode_solver_type")
  plot!(time_vs_dt,dts,times,xaxis=:log10,yaxis=:log10,label="$ode_solver_type")
  plot!(eu_vs_mem,memory,eus,xaxis=:log10,yaxis=:log10,label="$ode_solver_type")
end

# Finalize plots
plot!(eu_vs_dt,dts,1.0e-3 .* dts .^2,label="dt²",color=:black,ls=:dash)
plot!(eu_vs_dt,dts,1.0e-2 .* dts .^3,label="dt³",color=:black,ls=:dashdot)
eu_vs_dt_file = plotsdir("eu_vs_dt")
eu_vs_time_file = plotsdir("eu_vs_time")
time_vs_dt_file = plotsdir("time_vs_dt")
eu_vs_mem_file = plotsdir("eu_vs_memory")
savefig(eu_vs_dt,eu_vs_dt_file)
savefig(eu_vs_time,eu_vs_time_file)
savefig(eu_vs_mem,eu_vs_mem_file)
savefig(time_vs_dt,time_vs_dt_file)

end
