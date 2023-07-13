module MainDriver

using DrWatson
@quickactivate "GridapODEsTests.jl"

# Here you may include files from the source directory
include(srcdir("beltrami_flow.jl"))

println(
"""
Currently active project is: $(projectname())
Path of active project: $(projectdir())

"""
)

using .BeltramiFlow

params = BeltramiParams()
beltrami_flow_solver(params)

end
