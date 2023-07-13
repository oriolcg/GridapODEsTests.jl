module ManufacturedLinear
using Parameters
using Gridap
using DrWatson

export ManufacturedLinearParams
export manufactured_linear_solver

# """
# manufactured_linear_solver(params) → ||eᵤ||ₗ₂, ||eₚ||ₗ₂

# This function computes the L2-norm of the velocity and pressure error for
# the ManufacturedLinear flow problem.

# Manufactured lnear formulation:
# ```math
# u(x, y, t) = \begin{bmatrix}
# x\\
# −y
# \end{bmatrix} \sin\left(\pi/10 t\right)e^{t/25},
# p(x, y) = x + y.
# ```
# """
function manufactured_linear_solver(params)

  println("Executing test with the following settings:")
  println("-------------------------------------------")
  println(params)

  # Unpack variables
  @unpack vtk_output = params

  # Analytical solution
  u((x, y), t) = VectorValue(x,-y)*sin(π/10*t)*exp(t/25)
  u(t::Real) = x -> u(x,t)
  p((x, y)) = x+y

  # Discrete model
  @unpack ne = params
  𝒯 = CartesianDiscreteModel((0,1,0,1),(ne,ne))
  Ω = Interior(𝒯)
  Γ = Boundary(Ω)
  if vtk_output
    filename = datadir("sims","model")
    writevtk(𝒯,filename)
  end

  # FE spaces
  @unpack order = params
  refFEᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
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
  else
    error("ODE solver type not implemented")
  end

  # Solution
  xₕₜ = solve(ode_solver,op,xₕ₀,0.0,tf)

  # Postprocess
  global eᵤ,eₚ
  for (xₕ,t) in xₕₜ
    eᵤ = √( ∑( ∫( (u(t)-xₕ[1])⋅(u(t)-xₕ[1]) )dΩ ) )
    eₚ = √( ∑( ∫( (p-xₕ[2])⋅(p-xₕ[2]))dΩ ) )
    println("t = $t \n ========================")
  end

  println("Finished successfully!")
  return eᵤ,eₚ
end

"""
ManufacturedLinearParams

Struct with all required parameters in `manufactured_linear_solver`
"""
@with_kw struct ManufacturedLinearParams
  ode_solver_type::Symbol = :ThetaMethod
  ne::Integer = 2
  order::Integer = 2
  vtk_output::Bool = false
  ν::Real = 0.01
  dt::Real = 1.0e-2
  tf::Real = 1.0e-2
end

end
