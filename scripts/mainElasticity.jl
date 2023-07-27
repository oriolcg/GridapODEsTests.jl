module MainElasticityDriver

using DrWatson
using TimerOutputs
using Plots

@quickactivate "GridapODEsTests.jl"

# Here you may include files from the source directory
include(srcdir("beltrami_flow.jl"))
include(srcdir("manufactured_elasticity_linear.jl"))

println(
"""
Currently active project is: $(projectname())
Path of active project: $(projectdir())

"""
)

using .BeltramiFlow
using .ManufacturedElasticityLinear

# Initial run
params = ManufacturedElasticityLinearParams(vtk_output=true)
manufactured_elasticity_linear_solver(params)
params = ManufacturedElasticityLinearParams(ode_solver_type=:BackwardEuler)
println(manufactured_elasticity_linear_solver(params))
params = ManufacturedElasticityLinearParams(ode_solver_type=:GeneralizedAlpha)
manufactured_elasticity_linear_solver(params)
params = ManufacturedElasticityLinearParams(ode_solver_type=:RK_CN)
manufactured_elasticity_linear_solver(params)
params = ManufacturedElasticityLinearParams(ode_solver_type=:RK_BE)
println(manufactured_elasticity_linear_solver(params))
params = ManufacturedElasticityLinearParams(ode_solver_type=:RK_SDIRK)
manufactured_elasticity_linear_solver(params)
params = ManufacturedElasticityLinearParams(ode_solver_type=:RK_TRBDF2)
manufactured_elasticity_linear_solver(params)
params = ManufacturedElasticityLinearParams(ode_solver_type=:RK_ESDIRK3)
manufactured_elasticity_linear_solver(params)

# Timings
const to = TimerOutput()

# ODE solvers
ode_solver_types = [:ThetaMethod,:GeneralizedAlpha,:RK_CN,:RK_SDIRK,:RK_ESDIRK3]
tf = 0.5
dts = [tf/(2^i) for i in 0:3]

# Initialize plots
eu_vs_dt = plot(title="Error vs Δt",xlabel="Δt",ylabel="Error")
eu_vs_time = plot(title="Error vs CPU time",xlabel="CPU time",ylabel="Error")
time_vs_dt = plot(title="CPU time vs Δt",xlabel="Δt",ylabel="CPU time")
eu_vs_mem = plot(title="Error vs memory",xlabel="Memory",ylabel="Error")

# Loop over solvers
for ode_solver_type in ode_solver_types

  # Initialize output vectors
  eus = Float64[]
  times = Float64[]
  memory = Float64[]

  # Time convergence
  for dt in dts
    case = "$ode_solver_type"*"-$dt"
    params = ManufacturedElasticityLinearParams(tf=tf,dt=dt, ode_solver_type=ode_solver_type, ne=20, order=1)
    eᵤ, = @timeit to case manufactured_elasticity_linear_solver(params)
    push!(eus,eᵤ)
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
# plot!(eu_vs_dt,dts,1.0e-2 .* dts .^3,label="dt³",color=:black,ls=:dashdot)
eu_vs_dt_file = plotsdir("eu_vs_dt")
eu_vs_time_file = plotsdir("eu_vs_time")
time_vs_dt_file = plotsdir("time_vs_dt")
eu_vs_mem_file = plotsdir("eu_vs_memory")
savefig(eu_vs_dt,eu_vs_dt_file)
savefig(eu_vs_time,eu_vs_time_file)
savefig(eu_vs_mem,eu_vs_mem_file)
savefig(time_vs_dt,time_vs_dt_file)

end
