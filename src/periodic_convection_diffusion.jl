module PeriodicConvectionDiffusion
using Parameters
using Gridap
using DrWatson

export PeriodicConvectionDiffusionParams
export periodic_convection_diffusion_solver

"""
periodic_convection_diffusion_solver(params) → ||e||ₗ₂

This function computes the L2-norm of the error for the
PeriodicConvectionDiffusion problem:

```math
uₜ + (u⋅∇)u - νΔu = 0\qquad\text{in}\quad\Omega=[0,1]\times[0,1]
```
with periodic boundaries and initial condition:

```math
u(x, y, 0) = \begin{bmatrix}
\sin(2\pi (x+y))+0.005\cos(2\pi (64x + 63y)) \\
\sin(2\pi (x+y))+0.005\cos(2\pi (64x + 63y))
\end{bmatrix}.
```

The default parameters used in this test are based on the test in Section 5 of the paper
[1] Ascher, Uri M., Steven J. Ruuth, and Raymond J. Spiteri. "Implicit-explicit Runge-Kutta methods for time-dependent partial differential equations." Applied Numerical Mathematics 25.2-3 (1997): 151-167.
"""
function periodic_convection_diffusion_soler(params)

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
    writevtk(uₕ₀,filename)
  end

  # Measure
  degree = 2*order
  dΩ = Measure(Ω,degree)

  # Operators
  @unpack ν = params
  res(t,u,v) = ∫( (∂t(u) + (u⋅∇(u)))⋅v + ν*(∇(u)⊙∇(v)))dΩ
  op = TransientFEOperator(res,U,V)

  # Time stepping
  @unpack Δt = params
  t0 = 0.0
  tF = 0.5

  # ODE solver
  nls = NLSolver(show_trace=true,method=:newton,iterations=15)
  ode_solver = ThetaMethod(nls,Δt,θ=0.5)

  # Solution
  uₕₜ = solve(ode_solver,op,uₕ₀,t0,tF)

  # Postprocess
  e = 0.0
  for (uₕ,t) in uₕₜ
    e += √(∑(∫( (u(t)-uₕ)⋅(u(t)-uₕ) )dΩ))
    println("t = $t \n ========================")
  end

  # Return Output
  println("Finished successfully with error: $e")
  return e

end

"""
PeriodicConvectionDiffusionParams

This type is used to store the parameters of the PeriodicConvectionDiffusion
"""
@with_kw struct PeriodicConvectionDiffusionParams
  ν::Float64 = 0.01
  ne::Int64 = 128
  order::Int64 = 1
  Δt::Float64 = 0.00625
  vtk_output::Bool = false
end

end
