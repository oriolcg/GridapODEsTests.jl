module PeriodicConvectionDiffusion
using Parameters
using Gridap
using DrWatson

export PeriodicConvectionDiffusionParams
export periodic_convection_diffusion_solver

"""
periodic_convection_diffusion_solver(params) → ||u||ₗ₂

This function computes the L2-norm of the solution at the final time for the
PeriodicConvectionDiffusion problem:

```math
uₜ + (u⋅∇)u - νΔu = 0
```

with periodic boundaries and initial condition:

```math
u(x, y, 0) = [sin(2π (x+y)) + 0.005cos(2π (64x + 63y)),
              sin(2π (x+y)) + 0.005cos(2π (64x + 63y))].
```

The default parameters used in this test are based on the test in Section 5 of the paper
[1] Ascher, Uri M., Steven J. Ruuth, and Raymond J. Spiteri. "Implicit-explicit Runge-Kutta methods for time-dependent partial differential equations." Applied Numerical Mathematics 25.2-3 (1997): 151-167.

"""
function periodic_convection_diffusion_solver(params)

  println("Executing test with the following settings:")
  println("-------------------------------------------")
  println(params)

  # Unpack variables
  @unpack vtk_output = params

  # Analytical solution
  u((x, y), t) = VectorValue(sin(2π*(x+y))+0.005*cos(2π*(64x+63y)),
                             sin(2π*(x+y))+0.005*cos(2π*(64x+63y)))
  u(t::Real) = x -> u(x,t)

  # Discrete model
  @unpack ne = params
  𝒯 = CartesianDiscreteModel((0,1,0,1),(ne,ne),isperiodic=(true,true))
  Ω = Interior(𝒯)
  if vtk_output
    filename = datadir("sims","model_CD")
    writevtk(𝒯,filename)
  end

  # FE spaces
  @unpack order = params
  refFEᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
  V = TestFESpace(Ω,refFEᵤ,conformity=:H1)
  U = TransientTrialFESpace(V)

  # Initial solution
  uₕ₀ = interpolate_everywhere(u(0.0),U(0.0))

  # Output initial solution
  if vtk_output
    filename = datadir("sims","u0_CD")
    writevtk(Ω,filename,cellfields=["u0"=>uₕ₀])
  end

  # Measure
  degree = 2*order
  dΩ = Measure(Ω,degree)

  # Operators
  @unpack ν = params
  res(t,u,v) = ∫( (∂t(u) + (u⋅∇(u)))⋅v + ν*(∇(u)⊙∇(v)))dΩ
  lhs(t,u,v) = ∫( ∂t(u)⋅v )dΩ
  rhs(t,u,v) = ∫( -1.0*((u⋅∇(u)))⋅v - ν*(∇(u)⊙∇(v)))dΩ
  im_rhs(t,u,v) = ∫( - ν*(∇(u)⊙∇(v)))dΩ
  ex_rhs(t,u,v) = ∫( -1.0*((u⋅∇(u)))⋅v )dΩ
  op = TransientFEOperator(res,U,V)

  # Time stepping
  @unpack Δt,tf,ode_solver_type = params
  t0 = 0.0

  # ODE solver
  nls = NLSolver(show_trace=true,method=:newton,iterations=15)
  sol₀ = uₕ₀
  if ode_solver_type == :ThetaMethod
    ode_solver = ThetaMethod(nls,Δt,0.5)
  elseif ode_solver_type == :BackwardEuler
    ode_solver = BackwardEuler(nls,Δt)
  elseif ode_solver_type == :GeneralizedAlpha
    duₕ₀ = interpolate_everywhere(∂t(u)(0.0),U(0.0))
    sol₀ = (uₕ₀,duₕ₀)
    ode_solver = GeneralizedAlpha(nls,Δt,1.0)
  elseif ode_solver_type == :RK_BE
    op = TransientRungeKuttaFEOperator(lhs,rhs,U,V)
    ode_solver = RungeKutta(nls,LUSolver(),Δt,:BE_1_0_1)
  elseif ode_solver_type == :RK_CN
    op = TransientRungeKuttaFEOperator(lhs,rhs,U,V)
    ode_solver = RungeKutta(nls,LUSolver(),Δt,:CN_2_0_2)
  elseif ode_solver_type == :RK_SDIRK
    op = TransientRungeKuttaFEOperator(lhs,rhs,U,V)
    ode_solver = RungeKutta(nls,LUSolver(),Δt,:SDIRK_2_0_2)
  elseif ode_solver_type == :RK_ESDIRK3
    op = TransientRungeKuttaFEOperator(lhs,rhs,U,V)
    ode_solver = RungeKutta(nls,LUSolver(),Δt,:ESDIRK_3_1_2)
  elseif ode_solver_type == :IMEX_RK_FE_BE
    op = TransientIMEXRungeKuttaFEOperator(lhs,im_rhs,ex_rhs,U,V)
    ode_solver = IMEXRungeKutta(LUSolver(),LUSolver(),Δt,:IMEX_FE_BE_2_0_1)
  elseif ode_solver_type == :IMEX_RK_Midpoint
    op = TransientIMEXRungeKuttaFEOperator(lhs,im_rhs,ex_rhs,U,V)
    ode_solver = IMEXRungeKutta(LUSolver(),LUSolver(),Δt,:IMEX_Midpoint_2_0_2)
  else
    error("ODE solver type not implemented")
  end

  # Solution
  uₕₜ = solve(ode_solver,op,sol₀,t0,tf)

  # Postprocess
  global uₕ_final
  for (uₕ,t) in uₕₜ
    println("t = $t \n ========================")
    uₕ_final = √(∑(∫( uₕ⋅uₕ )dΩ))
  end

  return uₕ_final

end

"""
PeriodicConvectionDiffusionParams

This type is used to store the parameters of the PeriodicConvectionDiffusion
"""
@with_kw struct PeriodicConvectionDiffusionParams
  ode_solver_type::Symbol = :ThetaMethod
  ν::Float64 = 0.01
  ne::Int64 = 128
  order::Int64 = 1
  Δt::Float64 = 0.00625
  vtk_output::Bool = false
  tf::Real = 0.00625
end

end
