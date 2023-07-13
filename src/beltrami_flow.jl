module BeltramiFlow
using Parameters
using Gridap

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
  @unpack a, d = params

  # Analytical solution
  u((x, y, z), t) = VectorValue(
    −a * ( exp(ax) * sin(ay + dz) + exp(az) * cos(ax + dy) ),
    −a * ( exp(ay) * sin(az + dx) + exp(ax) * cos(ay + dz) ),
    −a * ( exp(az) * sin(ax + dy) + exp(ay) * cos(az + dx) )
  ) * exp(−d^2t)

  p((x, y, z)) = −a^2/2 * ( exp(2ax) + exp(2ay) + exp(2az) +
    2*sin(ax + dy) * cos(az + dx) * exp(a(y+z)) +
    2*sin(ay + dz) * cos(ax + dy) * exp(a(z+x)) +
    2*sin(az + dx) * cos(ay + dz) *exp(a(x+y)) )

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
end

end
