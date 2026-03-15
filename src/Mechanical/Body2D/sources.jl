
"""
    ForceAndTorque(; name, fx, fy, tau)

Force and torque applied to a frame.

# Connectors:

  - `frame` [Frame](@ref)
"""
@component function ForceAndTorque(; name)
    pars = @parameters begin
    end

    systems = @named begin
        fx = RealInput()
        fy = RealInput()
        tau = RealInput()
        frame = Frame()
    end

    vars = @variables begin
    end

    equations = Equation[
        frame.fx ~ -fx.u, # TODO: Why would forces and torque be applied negative (also used negative in rotational torque source, but why?)
        frame.fy ~ -fy.u,
        frame.tau ~ -tau.u,
    ]

    return System(equations, t, [], pars; name, systems)
end

@component function SupportedTorque(; name)
    pars = @parameters begin
    end

    systems = @named begin
        tau = RealInput()
        frame = Frame()
        support = Frame()
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

    return System(equations, t, [], pars; name, systems)
end
