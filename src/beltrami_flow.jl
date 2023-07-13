module BeltramiFlow
using Parameters
using Gridap
using DrWatson

export BeltramiParams
export beltrami_flow_solver

# """
# beltrami_flow_solver(params) → ||eᵤ||ₗ₂, ||eₚ||ₗ₂

# This function computes the L2-norm of the velocity and pressure error for
# the Beltrami flow problem.

# Beltrami flow formulation:
# ```math
# u(x, y, z, t) = \begin{bmatrix}
# −a \left[e^{ax} \sin(ay + dz) + e^{az} \cos(ax + dy)\right]\\
# −a \left[e^{ay} \sin(az + dx) + e^{ax} \cos(ay + dz)\right]\\
# −a \left[e^{az} \sin(ax + dy) + e^{ay} \cos(az + dx)\right]
# \end{bmatrix} e^{−d^2t},
# p(x, y, z) = −\frac{a^2}{2} \left[ e^{2ax} + e^{2ay} + e^{2az} +
#   2\sin(ax + dy) \cos(az + dx) e^{a(y+z)} +
#   2\sin(ay + dz) \cos(ax + dy) e^{a(z+x)} +
#   2\sin(az + dx) \cos(ay + dz) e^{a(x+y)}.
# ```
# """
function beltrami_flow_solver(params)

  # Unpack variables
  @unpack a, d, vtk_output = params

  # Analytical solution
  u((x, y, z), t) = VectorValue(
    −a * ( exp(a*x) * sin(a*y + d*z) + exp(a*z) * cos(a*x + d*y) ),
    −a * ( exp(a*y) * sin(a*z + d*x) + exp(a*x) * cos(a*y + d*z) ),
    −a * ( exp(a*z) * sin(a*x + d*y) + exp(a*y) * cos(a*z + d*x) )
  ) * exp(−d^2t)
  u(t::Real) = x -> u(x,t)
  p((x, y, z)) = −a^2/2 * ( exp(2a*x) + exp(2a*y) + exp(2a*z) +
    2*sin(a*x + d*y) * cos(a*z + d*x) * exp(a*(y+z)) +
    2*sin(a*y + d*z) * cos(a*x + d*y) * exp(a*(z+x)) +
    2*sin(a*z + d*x) * cos(a*y + d*z) * exp(a*(x+y)) )

  # Discrete model
  @unpack ne = params
  𝒯 = CartesianDiscreteModel((-1,1,-1,1,-1,1),(ne,ne,ne))
  Ω = Interior(𝒯)
  Γ = Boundary(Ω)
  if vtk_output
    filename = datadir("sims","model")
    writevtk(𝒯,filename)
  end

  # FE spaces
  @unpack order = params
  refFEᵤ = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
  refFEₚ = ReferenceFE(lagrangian,Float64,order-1)
  V = TestFESpace(Ω,refFEᵤ,conformity=:H1,dirichlet_tags="boundary")
  Q = TestFESpace(Ω,refFEₚ,conformity=:C0,constraint=:zeromean)
  U = TransientTrialFESpace(V,u)
  P = TrialFESpace(Q)
  Y = MultiFieldFESpace([V,Q])
  X = TransientMultiFieldFESpace([U,P])

  # Initial solution
  xₕ₀ = interpolate_everywhere([u(0.0),p],X(0.0))

  # Output Initial solution
  if vtk_output
    filename = datadir("sims","xh0")
    writevtk(Ω,filename,cellfields=["u"=>xₕ₀[1],"p"=>xₕ₀[2]])
  end

  # Measures
  dΩ = Measure(Ω,2*order)

  # Weak form
  @unpack ν = params
  b(a,u,v) = 0.5*( (a⋅∇(u))⋅v - (a⋅∇(v))⋅u )
  # res(t,(u,p),(v,q)) = ∫( ∂t(u)⋅v + ν*(∇(u)⊙∇(v)) - p*(∇⋅v) + q*(∇⋅u) )dΩ
  res(t,(u,p),(v,q)) = ∫( ∂t(u)⋅v + b(u,u,v) + ν*(∇(u)⊙∇(v)) - p*(∇⋅v) + q*(∇⋅u) )dΩ
  op = TransientFEOperator(res,X,Y)

  # ODE solver
  @unpack ode_solver_type, dt, tf = params
  nls = NLSolver(show_trace=true,method=:newton,iterations=15)
  if ode_solver_type == :ThetaMethod
    ode_solver = ThetaMethod(nls,dt,0.5)
  end

  # Solution
  xₕₜ = solve(ode_solver,op,xₕ₀,0.0,tf)

  # Postprocess
  for (xₕ,t) in xₕₜ
    eᵤ = √( ∑( ∫( (u(t)-xₕ[1])⋅(u(t)-xₕ[1]) )dΩ ) )
    eₚ = √( ∑( ∫( (p-xₕ[2])⋅(p-xₕ[2]) )dΩ ) )
    println("t = $t\n ========================")
    println("eᵤ = $eᵤ")
    println("eₚ = $eₚ")
  end

  println("Finished successfully!")
  return nothing
end

"""
BeltramiParams

Struct with all required parameters in `beltrami_flow_solver`
"""
@with_kw struct BeltramiParams
  ode_solver_type::Symbol = :ThetaMethod
  a::Real = π/4
  d::Real = π/2
  ne::Integer = 2
  order::Integer = 2
  vtk_output::Bool = false
  ν::Real = 0.01
  dt::Real = 1.0e-2
  tf::Real = 2.0e-2
end

end
