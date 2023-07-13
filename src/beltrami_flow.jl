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
  Q = TestFESpace(Ω,refFEₚ,conformity=:L2,constraint=:zeromean)
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

  println("Finished successfully!")
  return nothing
end

"""
BeltramiParams

Struct with all required parameters in `beltrami_flow_solver`
"""
@with_kw struct BeltramiParams
  ode_solver::Symbol = :ThetaMethod
  a::Real = π/4
  d::Real = π/2
  ne::Integer = 2
  order::Integer = 2
  vtk_output::Bool = false
end

end
