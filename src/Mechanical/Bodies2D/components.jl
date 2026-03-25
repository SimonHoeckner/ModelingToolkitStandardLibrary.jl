
@component function Fixed2d(; name, x, y, phi, k = 1e6)
    pars = @parameters begin
        x = x
        y = y
        phi = phi
    end

    vars = @variables begin
        rx(t), [guess = 0.0]
        ry(t), [guess = 0.0]
    end

    systems = @named begin
        frame = Frame2d()
    end

    equations = [
        frame.x ~ rx,
        frame.y ~ ry,
        frame.phi ~ phi,
        frame.fx ~ -k * (rx - x),
        frame.fy ~ -k * (ry - y),
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Rigidbody2d(;name, m, J, x, y, phi, vx, vy, w)

Body in 2 dimensions with a mass and a moment of inertia. By itself, it can move and rotate freely.

# Connectors:

  - `frame` [Frame2d](@ref)

# Parameters:

  - `m`: [`kg`] Mass of the body
  - `J`: [`kg·m²`] Moment of inertia of the body
  - `g`: [`m/s²`] Gravitational acceleration of the body in negative y-direction

# States:

  - `x(t)`: [`m`] Absolute x-position of the body's center of mass
  - `y(t)`: [`m`] Absolute y-position of the body's center of mass
  - `phi(t)`: [`rad`] Absolute rotation angle of the body
  - `vx(t)`: [`m/s`] Absolute x-velocity of the body's center of mass
  - `vy(t)`: [`m/s`] Absolute y-velocity of the body's center of mass
  - `w(t)`: [`rad/s`] Absolute angular velocity of body
  - `ax(t)`: [`m/s²`] Absolute x-acceleration of the body's center of mass
  - `ay(t)`: [`m/s²`] Absolute y-acceleration of the body's center of mass
  - `a(t)`: [`rad/s²`] Absolute angular acceleration of the body
"""
@component function Rigidbody2d(; m = nothing, J = nothing, g = 9.81,
                                  x = nothing, y = nothing, phi = nothing,
                                  vx = nothing, vy = nothing, w = nothing,
                                  name)
    @symcheck m > 0 || throw(ArgumentError("Expected `m` to be positive"))
    @symcheck J > 0 || throw(ArgumentError("Expected `J` to be positive"))

    pars = @parameters begin
        m = m, [description = "Mass of the body"]
        J = J, [description = "Moment of inertia of the body"]
        g = g, [description = "Gravitational acceleration acting in negative y-direction"]
    end

    systems = @named begin
        frame = Frame2d()
    end

    vars = @variables begin
        x(t) = x, [description = "Absolute x-position of the body's center of mass", guess = 0.0]
        y(t) = y, [description = "Absolute y-position of the body's center of mass", guess = 0.0]
        phi(t) = phi, [description = "Absolute rotation angle of the body", guess = 0.0]
        vx(t) = vx, [description = "Absolute x-velocity of the body's center of mass", guess = 0.0]
        vy(t) = vy, [description = "Absolute y-velocity of the body's center of mass", guess = 0.0]
        w(t) = w, [description = "Absolute angular velocity of body", guess = 0.0]
        ax(t), [description = "Absolute x-acceleration of the body's center of mass", guess = 0.0]
        ay(t), [description = "Absolute y-acceleration of the body's center of mass", guess = 0.0]
        a(t), [description = "Absolute angular acceleration of the body", guess = 0.0]
    end

    equations = Equation[
        frame.x ~ x,
        frame.y ~ y,
        frame.phi ~ phi,
        D(x) ~ vx,
        D(y) ~ vy,
        D(phi) ~ w,
        D(vx) ~ ax,
        D(vy) ~ ay,
        D(w) ~ a,
        m * ax ~ frame.fx,
        m * ay ~ frame.fy - m * g,
        J * a ~ frame.tau,
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Frame2dOffset(; name, x_0, y_0, k)

Rigid link between `frame_a` and `frame_b` with offset `x_0` and `y_0` in the local coordinate system of `frame_a`.
The offset is implemented as a spring with high stiffness `k` to avoid numerical issues.

# Connectors:

  - `frame_a` [Frame2d](@ref)
  - `frame_b` [Frame2d](@ref)

# Parameters:

  - `x_0`: [`m`] x-offset of `frame_b` relative to `frame_a` in the local coordinates of `frame_a`
  - `y_0`: [`m`] y-offset of `frame_b` relative to `frame_a` in the local coordinates of `frame_a`
  - `k`: [`N/m`] Spring stiffness
"""
@component function Frame2dOffset(; name, x_0, y_0, k = 1e6)
    pars = @parameters begin
    end

    systems = @named begin
        frame_a = Frame2d()
        frame_b = Frame2d()
    end

    rx_0 = cos(frame_a.phi) * x_0 - sin(frame_a.phi) * y_0
    ry_0 = sin(frame_a.phi) * x_0 + cos(frame_a.phi) * y_0

    vars = @variables begin
        rx(t) = rx_0
        ry(t) = ry_0
        λx(t), [guess = 0.0]
        λy(t), [guess = 0.0]
    end

    equations = [
        # Geometry:
        frame_b.x ~ frame_a.x + rx,
        frame_b.y ~ frame_a.y + ry,
        frame_b.phi ~ frame_a.phi,
        # Torque balance:
        frame_a.tau + frame_b.tau + rx * frame_b.fy - ry * frame_b.fx ~ 0,
        # Spring forces:
        λx ~ k * (rx - rx_0),
        λy ~ k * (ry - ry_0),
        frame_a.fx ~ -λx,
        frame_b.fx ~ λx,
        frame_a.fy ~ -λy,
        frame_b.fy ~ λy,
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    RevoluteJoint2d(; name, phi, k)

A joint between two [Frame2d](@ref) with free rotation between the two frames.

# Connectors:

  - `frame_a` [Frame2d](@ref)
  - `frame_b` [Frame2d](@ref)
"""
@component function RevoluteJoint2d(; name, phi = nothing, k = 1e6)
    pars = @parameters begin
    end

    vars = @variables begin
        phi(t) = phi
        rx(t), [guess = 0.0]
        ry(t), [guess = 0.0]
        λx(t), [guess = 0.0]
        λy(t), [guess = 0.0]
    end

    systems = @named begin
        frame_a = Frame2d()
        frame_b = Frame2d()
    end

    equations = [
        frame_b.x ~ frame_a.x + rx,
        frame_b.y ~ frame_a.y + ry,
        frame_a.fx + frame_b.fx ~ 0.0,
        frame_a.fy + frame_b.fy ~ 0.0,
        λx ~ k * rx,
        λy ~ k * ry,
        frame_a.fx ~ -λx,
        frame_b.fx ~ λx,
        frame_a.fy ~ -λy,
        frame_b.fy ~ λy,

        frame_b.phi - frame_a.phi ~ phi,
        frame_a.tau ~ 0.0,
        frame_b.tau ~ 0.0,
    ]

    return System(equations, t, vars, pars; name, systems)
end

@component function PrismaticJoint2d(; name, l = nothing, phi, k = 1e6)
    pars = @parameters begin
    end

    vars = @variables begin
        l(t) = l
        epsilon(t) = 0.0
    end

    systems = @named begin
        frame_a = Frame2d()
        frame_b = Frame2d()
    end

    joint_angle = phi + frame_a.phi

    f_tangential_a = -frame_a.fx * sin(joint_angle) + frame_a.fy * cos(joint_angle)
    f_tangential_b = -frame_b.fx * sin(joint_angle) + frame_b.fy * cos(joint_angle)

    eqs = [
        # Geometry:
        frame_a.phi ~ frame_b.phi,
        (frame_b.x - frame_a.x) * sin(joint_angle) - (frame_b.y - frame_a.y) * cos(joint_angle) ~ epsilon,
        (frame_b.x - frame_a.x) * cos(joint_angle) + (frame_b.y - frame_a.y) * sin(joint_angle) ~ l,

        # Forces:
        f_tangential_a ~ k * epsilon,  # Forces tangential to the joint are transmitted
        f_tangential_b ~ -k * epsilon, # Forces tangential to the joint are transmitted
        frame_a.fx * cos(joint_angle) + frame_a.fy * sin(joint_angle) ~ 0.0, # Forces in direction of the joint are not transmitted
        frame_b.fx * cos(joint_angle) + frame_b.fy * sin(joint_angle) ~ 0.0, # Forces in direction of the joint are not transmitted
        # Torque is transmitted
        frame_a.tau + frame_b.tau + f_tangential_b * l ~ 0.0,
    ]

    return System(eqs, t, vars, pars; name, systems)
end
