"""
Library to model 2-dimensional mechanical systems with rigid bodies and links of different kinds between them.

This is very non-functioning and will probably never be.
"""
module Body2D

using ModelingToolkitBase, Symbolics, IfElse
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using ...Blocks: RealInput, RealOutput
import ...@symcheck

export Frame
include("utils.jl")

export Body2d, FrameOffset
include("components.jl")

export ForceAndTorque, SupportedTorque
include("sources.jl")

end
