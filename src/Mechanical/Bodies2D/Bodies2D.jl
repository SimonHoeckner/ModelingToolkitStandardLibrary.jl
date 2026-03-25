"""
Library to model 2-dimensional mechanical systems with rigid bodies and links of different kinds between them.
"""
module Bodies2D

using ModelingToolkitBase, Symbolics, IfElse
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using ...Blocks: RealInput, RealOutput
import ...@symcheck

export Frame2d
include("utils.jl")

export Fixed2d, Rigidbody2d, Frame2dOffset, RevoluteJoint2d, PrismaticJoint2d
include("components.jl")

export ForceAndTorque2d, SupportedTorque2d
include("sources.jl")

end
