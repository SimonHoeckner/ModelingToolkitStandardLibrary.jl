
"""
    Frame2d(;name, x, y, phi)

Frame of reference and attachment point for 2D rigid body mechanics.

# States:

  - `x(t)`: [`m`] Absolute x-position of the frame
  - `y(t)`: [`m`] Absolute y-position of the frame
  - `phi(t)`: [`rad`] Absolute rotation angle of the frame
  - `fx(t)`: [`N`] x-component of the applied force
  - `fy(t)`: [`N`] y-component of the applied force
  - `tau(t)`: [`N·m`] Applied torque
"""
@connector function Frame2d(; name, x = nothing, y = nothing, phi = nothing, fx = nothing, fy = nothing, tau = nothing)
    vars = @variables begin
        x(t) = x, [description = "Absolute x-position of the frame", guess = 0.0]
        y(t) = y, [description = "Absolute y-position of the frame", guess = 0.0]
        phi(t) = phi, [description = "Absolute rotation angle of the frame", guess = 0.0]
        fx(t), [description = "x-component of the applied force", connect = Flow, guess = 0.0]
        fy(t), [description = "y-component of the applied force", connect = Flow, guess = 0.0]
        tau(t), [description = "Applied torque", connect = Flow, guess = 0.0]
    end

    return System(Equation[], t, vars, []; name)
end
