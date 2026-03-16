
"""
    Rigidbody2d(;name, m, J, x, y, phi, vx, vy, w)

Body in 2 dimensions with a mass and a moment of inertia. By itself, it can move and rotate freely.

# Connectors:

  - `frame` [Frame2d](@ref)

# Parameters:

  - `m`: [`kg`] Mass of the body
  - `J`: [`kg·m²`] Moment of inertia of the body

# States:

  - `x(t)`: [`m`] Absolute x-position of the body's center of mass
  - `y(t)`: [`m`] Absolute y-position of the body's center of mass
  - `phi(t)`: [`rad`] Absolute rotation angle of the body
  - `vx(t)`: [`m/s`] Absolute x-velocity of the body's center of mass
  - `vy(t)`: [`m/s`] Absolute y-velocity of the body's center of mass
  - `w(t)`: [`rad/s`] Absolute angular velocity of body
"""
@component function Rigidbody2d(; m = nothing, J = nothing,
                             x = nothing, y = nothing, phi = nothing,
                             vx = nothing, vy = nothing, w = nothing,
                             name)
    @symcheck m > 0 || throw(ArgumentError("Expected `m` to be positive"))
    @symcheck J > 0 || throw(ArgumentError("Expected `J` to be positive"))

    pars = @parameters begin
        m = m, [description = "Mass of the body"]
        J = J, [description = "Moment of inertia of the body"]
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
    end

    equations = Equation[
        frame.x ~ x,
        frame.y ~ y,
        frame.phi ~ phi,
        D(x) ~ vx,
        D(y) ~ vy,
        D(phi) ~ w,
        m * D(vx) ~ frame.fx,
        m * D(vy) ~ frame.fy,
        J * D(w) ~ frame.tau,
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
        x_0 = x_0
        y_0 = y_0
    end

    vars = @variables begin
        rx(t) = x_0
        ry(t) = y_0
        λx(t), [guess = 0.0]
        λy(t), [guess = 0.0]
    end

    systems = @named begin
        frame_a = Frame2d()
        frame_b = Frame2d()
    end

    rx_0 = cos(frame_a.phi) * x_0 - sin(frame_a.phi) * y_0
    ry_0 = sin(frame_a.phi) * x_0 + cos(frame_a.phi) * y_0

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
