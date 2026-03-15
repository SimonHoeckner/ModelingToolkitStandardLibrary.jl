
"""
    Body2d(;name, m, J, x, y, phi, vx, vy, w)

Body in 2 dimensions with a mass and a moment of inertia. By itself, it can move and rotate freely.

# Connectors:

  - `frame` [Frame](@ref)

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
  - `ax(t)`: [`m/s²`] Absolute x-acceleration of the body's center of mass
  - `ay(t)`: [`m/s²`] Absolute y-acceleration of the body's center of mass
  - `a(t)`: [`rad/s²`] Absolute angular acceleration of the body
"""
@component function Body2d(; m = nothing, J = nothing,
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
        frame = Frame()
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
        m * ay ~ frame.fy,
        J * a ~ frame.tau,
    ]

    return System(equations, t, vars, pars; name, systems)
end

@component function FrameOffset(; name, x_0, y_0, k = 1e6)
    L_0 = sqrt(x_0^2 + y_0^2)

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
        frame_a = Frame()
        frame_b = Frame()
    end

    rx_0 = cos(frame_a.phi) * x_0 - sin(frame_a.phi) * y_0
    ry_0 = sin(frame_a.phi) * x_0 + cos(frame_a.phi) * y_0

    equations = [
        # Geometry:
        frame_b.x ~ frame_a.x + rx,
        frame_b.y ~ frame_a.y + ry,
        frame_b.phi ~ frame_a.phi,
        # Force balance:
        # frame_a.fx + frame_b.fx ~ λx * rx / rx_0,
        # frame_a.fy + frame_b.fy ~ λy * ry / ry_0,
        frame_a.fx ~ -λx,
        frame_b.fx ~ λx,
        frame_a.fy ~ -λy,
        frame_b.fy ~ λy,
        # Torque balance:
        frame_a.tau + frame_b.tau + rx * frame_b.fy - ry * frame_b.fx ~ 0,
        # Compliance:
        λx ~ k * (rx - rx_0),
        λy ~ k * (ry - ry_0),
    ]

    return System(equations, t, vars, pars; name, systems)
end
