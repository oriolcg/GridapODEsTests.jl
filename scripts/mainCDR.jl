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
uₕ₀ = periodic_convection_diffusion_solver(params)
params = PeriodicConvectionDiffusionParams(ne=10,Δt=0.0001,tf = 0.01)
uₕ_ref = periodic_convection_diffusion_solver(params)

# Timings
const to = TimerOutput()

# ODE solvers
ode_solver_types = [:ThetaMethod,:GeneralizedAlpha,:RK_CN,:RK_SDIRK,:RK_ESDIRK3,:IMEX_RK_FE_BE,:IMEX_RK_Midpoint]
tf = 0.01
dts = [tf/(2^i) for i in 0:3]

# Initialize plots
eu_vs_dt = plot(title="Error vs Δt",xlabel="Δt",ylabel="Error")
eu_vs_time = plot(title="Error vs CPU time",xlabel="CPU time",ylabel="Error")
time_vs_dt = plot(title="CPU time vs Δt",xlabel="Δt",ylabel="CPU time")
eu_vs_mem = plot(title="Error vs memory",xlabel="Memory",ylabel="Error")

# Loop over solvers
for ode_solver_type in ode_solver_types

  # Initialize output vectors
  es = Float64[]
  times = Float64[]
  memory = Float64[]

  # Time convergence
  for dt in dts
    case = "$ode_solver_type"*"-$dt"
    params = PeriodicConvectionDiffusionParams(tf=tf,Δt=dt, ode_solver_type=ode_solver_type, ne=10, order=1)
    uₕ = @timeit to case periodic_convection_diffusion_solver(params)
    e = abs(uₕ-uₕ_ref)
    push!(es,e)
    push!(times,TimerOutputs.time(to[case]))
    push!(memory,TimerOutputs.allocated(to[case]))
  end

  # Plotting
  plot!(eu_vs_dt,dts,es,xaxis=:log10,yaxis=:log10,label="$ode_solver_type")
  plot!(eu_vs_time,times,es,xaxis=:log10,yaxis=:log10,label="$ode_solver_type")
  plot!(time_vs_dt,dts,times,xaxis=:log10,yaxis=:log10,label="$ode_solver_type")
  plot!(eu_vs_mem,memory,es,xaxis=:log10,yaxis=:log10,label="$ode_solver_type")

  println(es)
end

# Finalize plots
plot!(eu_vs_dt,dts,1.0e-3 .* dts .^2,label="dt²",color=:black,ls=:dash)
# plot!(eu_vs_dt,dts,1.0e-2 .* dts .^3,label="dt³",color=:black,ls=:dashdot)
eu_vs_dt_file = plotsdir("CDR","eu_vs_dt")
eu_vs_time_file = plotsdir("CDR","eu_vs_time")
time_vs_dt_file = plotsdir("CDR","time_vs_dt")
eu_vs_mem_file = plotsdir("CDR","eu_vs_memory")
savefig(eu_vs_dt,eu_vs_dt_file)
savefig(eu_vs_time,eu_vs_time_file)
savefig(eu_vs_mem,eu_vs_mem_file)
savefig(time_vs_dt,time_vs_dt_file)

end
