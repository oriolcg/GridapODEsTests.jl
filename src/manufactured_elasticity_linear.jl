module ManufacturedElasticityLinear
using Parameters
using Gridap
using DrWatson

export ManufacturedElasticityLinearParams
export manufactured_elasticity_linear_solver

# """
# manufactured_elasticity_linear_solver(params) → ||eᵤ||ₗ₂, ||eₚ||ₗ₂

# This function computes the L2-norm of the displacement error for
# the ManufacturedElasticityLinear problem.

# Manufactured linear formulation in (velocity as unknown):
# ```math
# u(x, y, t) = \begin{bmatrix}
# x\\
# −y
# \end{bmatrix} \sin\left(\pi/10 t\right)e^{t/25}.
# ```
# """
function manufactured_elasticity_linear_solver(params)

  println("Executing test with the following settings:")
  println("-------------------------------------------")
  println(params)

  # Unpack variables
  @unpack vtk_output = params
  @unpack ρ,ν,E = params
  μ = E/(2(1+ν))
  λ = E*ν/((1+ν)*(1-2ν))

  # Analytical solution
  I = TensorValue(1.0,0.0,0.0,1.0)
  ux((x, y))= VectorValue(x,-y)
  σux(x) = μ * (∇(ux)(x) + ∇(ux)(x)' ) + λ * (∇⋅ux)(x) * I
  ut(t)::Real = sin(π/10*t)*exp(t/25)
  dut(t)::Real = π/10*cos(π/10*t)exp(t/25) + 1/25*sin(π/10*t)*exp(t/25)
  u(x, t) = ux(x)*ut(t)
  u(t::Real) = x -> u(x,t)

  # Discrete model
  @unpack ne = params
  𝒯 = CartesianDiscreteModel((0,1,0,1),(ne,ne))
  Ω = Interior(𝒯)
  Γ = Boundary(Ω)
  n = get_normal_vector(Γ)
  if vtk_output
    filename = datadir("sims","model")
    writevtk(𝒯,filename)
  end

  # FE spaces
  @unpack order = params
  refFEᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
  V = TestFESpace(Ω,refFEᵤ,conformity=:H1)#,dirichlet_tags="boundary")
  U = TransientTrialFESpace(V)#,u)

  # Initial solution
  uₕ₀ = interpolate_everywhere(u(0.0),U(0.0))

  # Source terms
  f(t) = x -> dut(t)*ux(x) - 1/ρ*(∇⋅σux)(x)*ut(t)
  g(t) = x -> (σux(x)⋅n)*ut(t)

  # Output Initial solution
  if vtk_output
    filename = datadir("sims","manufacturedElasticityLinear_xh0")
    writevtk(Ω,filename,cellfields=["u"=>uₕ₀])
  end

  # Measures
  dΩ = Measure(Ω,2*order)
  dΓ = Measure(Γ,2*order)

  # Weak form
  σ(u) = 2μ*ε(u) + λ*(∇⋅u)*I
  res(t,u,v) = ∫( ρ*∂t(u)⋅v + σ(u)⊙ε(v) - ρ*(f(t)⋅v) )dΩ - ∫(σux⊙(n⊗v)*ut(t))dΓ
  lhs(t,u,v) = ∫( ρ*∂t(u)⋅v )dΩ
  rhs(t,u,v) = ∫( - σ(u)⊙ε(v) + ρ*(f(t)⋅v) )dΩ + ∫(σux⊙(n⊗v)*ut(t))dΓ
  op = TransientFEOperator(res,U,V)

  # ODE solver
  @unpack ode_solver_type, dt, tf = params
  nls = NLSolver(show_trace=true,method=:newton,iterations=15)
  sol₀ = uₕ₀
  if ode_solver_type == :ThetaMethod
    ode_solver = ThetaMethod(nls,dt,0.5)
  elseif ode_solver_type == :BackwardEuler
      ode_solver = BackwardEuler(nls,dt)
  elseif ode_solver_type == :GeneralizedAlpha
    duₕ₀ = interpolate_everywhere(∂t(u)(0.0),U(0.0))
    sol₀ = (uₕ₀,duₕ₀)
    ode_solver = GeneralizedAlpha(nls,dt,1.0)
  elseif ode_solver_type == :RK_BE
    op = TransientRungeKuttaFEOperator(lhs,rhs,U,V)
    ode_solver = RungeKutta(nls,LUSolver(),dt,:BE_1_0_1)
  elseif ode_solver_type == :RK_CN
    op = TransientRungeKuttaFEOperator(lhs,rhs,U,V)
    ode_solver = RungeKutta(nls,LUSolver(),dt,:CN_2_0_2)
  elseif ode_solver_type == :RK_SDIRK
    op = TransientRungeKuttaFEOperator(lhs,rhs,U,V)
    ode_solver = RungeKutta(nls,LUSolver(),dt,:SDIRK_2_0_2)
  elseif ode_solver_type == :RK_TRBDF2
    op = TransientRungeKuttaFEOperator(lhs,rhs,U,V)
    ode_solver = RungeKutta(nls,LUSolver(),dt,:TRBDF2_3_2_3)
  elseif ode_solver_type == :RK_ESDIRK3
    op = TransientRungeKuttaFEOperator(lhs,rhs,U,V)
    ode_solver = RungeKutta(nls,LUSolver(),dt,:ESDIRK_3_1_2)
  else
    error("ODE solver type not implemented")
  end

  # Solution
  uₕₜ = solve(ode_solver,op,sol₀,0.0,tf)

  # Postprocess
  eᵤ=Float64[]
  eₚ=Float64[]
  for (uₕ,t) in uₕₜ
    push!(eᵤ,√( ∑( ∫( (u(t)-uₕ)⋅(u(t)-uₕ) )dΩ ) ))
    println("t = $t \n ========================")
  end

  println("Finished successfully!")
  return (last(eᵤ),)
end

"""
ManufacturedElasticityLinearParams

Struct with all required parameters in `manufactured_elasticity_linear_solver`
"""
@with_kw struct ManufacturedElasticityLinearParams
  ode_solver_type::Symbol = :ThetaMethod
  ne::Integer = 2
  order::Integer = 2
  vtk_output::Bool = false
  ρ::Real = 0.5
  E::Real = 1e4
  ν::Real = 0.3
  dt::Real = 1.0e-2
  tf::Real = 1.0e-2
end

end
