"""
Library of mechanical models.
"""
module Mechanical

using ModelingToolkitBase

include("Rotational/Rotational.jl")
include("Translational/Translational.jl")
include("TranslationalPosition/TranslationalPosition.jl")
include("TranslationalModelica/TranslationalModelica.jl")
include("MultiBody2D/MultiBody2D.jl")
include("Bodies2D/Bodies2D.jl")

end
