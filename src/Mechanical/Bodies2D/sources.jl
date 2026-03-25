
"""
    ForceAndTorque2d(; name)

The provided force and torque are applied to a frame.

# Connectors:

  - `fx` [RealInput](@ref)
  - `fy` [RealInput](@ref)
  - `tau` [RealInput](@ref)
  - `frame` [Frame2d](@ref)
"""
@component function ForceAndTorque2d(; name)
    pars = @parameters begin
    end

    systems = @named begin
        fx = RealInput()
        fy = RealInput()
        tau = RealInput()
        frame = Frame2d()
    end

    vars = @variables begin
    end

    equations = Equation[
        frame.fx ~ -fx.u, # TODO: Why would forces and torque be applied negative (also used negative in rotational torque source, but why?)
        frame.fy ~ -fy.u,
        frame.tau ~ -tau.u,
    ]

    return System(equations, t, vars, pars; name, systems)
end


"""
    SupportedTorque2d(; name)

Torque applied to a frame. The same torque is applied to the support frame in opposite direction.

# Connectors:

  - `tau` [RealInput](@ref)
  - `frame` [Frame2d](@ref)
  - `support` [Frame2d](@ref)
"""
@component function SupportedTorque2d(; name)
    pars = @parameters begin
    end

    systems = @named begin
        tau = RealInput()
        frame = Frame2d()
        support = Frame2d()
    end

    vars = @variables begin
    end

    equations = Equation[
        frame.tau ~ -tau.u, # TODO: Why exactly these signs?
        support.tau ~ tau.u,
        frame.fx ~ 0.0,
        frame.fy ~ 0.0,
        support.fx ~ 0.0,
        support.fy ~ 0.0
    ]

    return System(equations, t, vars, pars; name, systems)
end
