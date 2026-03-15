using ModelingToolkitStandardLibrary.Mechanical.Body2D
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit, OrdinaryDiffEq, SciMLBase, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "free body, constant velocity" begin
    # A body with an initial x-velocity of 2m/s and no external forces.
    # Expected: x-position of 2m/s * 3s after 3 seconds of simulation.

    @named sys = Body2d(m = 1.0, J = 1.0, x = 0.0, y = 0.0, phi = 0.0, vx = 2.0, vy = 0.0, w = 0.0)
    sys = mtkcompile(sys)
    prob = ODEProblem(sys, [], (0.0, 5.0))
    sol = solve(prob, Tsit5())

    @test SciMLBase.successful_retcode(sol)
    @test sol[sys.x][end] ≈ 10.0
    @test sol[sys.y][end] ≈ 0.0
    @test sol[sys.phi][end] ≈ 0.0
    @test all(sol[sys.vx] .≈ 2.0)
end

@testset "free body, constant force" begin
    # A body with mass 2 kg is subjected to constant a force in x-direction of 4 N.
    # Expected: x-position of 1/2 * (4N / 2kg) * 3s = 3m after 3 seconds of simulation.
    function ForcedBody(;name)
        @named body = Body2d(m = 2.0, J = 1.0, x = 0.0, y = 0.0, phi = 0.0,
                                vx = 0.0, vy = 0.0, w = 0.0)
        @named force = ForceAndTorque()
        @named force_step = Blocks.Step(height = 4.0, start_time = 0.0)
        @named zero_fy = Blocks.Step(height = 0.0, start_time = 0.0)
        @named zero_tau = Blocks.Step(height = 0.0, start_time = 0.0)
        eqs = [
            connect(body.frame, force.frame),
            connect(force.fx, force_step.output),
            connect(force.fy, zero_fy.output),
            connect(force.tau, zero_tau.output)
        ]
        return System(eqs, t, [], []; name, systems = [body, force, force_step, zero_fy, zero_tau])
    end

    @named sys = ForcedBody()
    sys = mtkcompile(sys)
    prob = ODEProblem(sys, [], (0.0, 3.0))
    sol = solve(prob, Tsit5())

    @test SciMLBase.successful_retcode(sol)
    @test sol[sys.body.x][end] ≈ 0.5 * 2.0 * 3.0^2 atol = 1e-3
    @test sol[sys.body.vx][end] ≈ 2.0 * 3.0 atol = 1e-3
end

@testset "two connected 2d bodies" begin
    m1 = 1.0 # kg
    m2 = 0.2 # kg
    J1 = 5.0 # kg·m²
    J2 = 2.0 # kg·m²
    l = 1.0 # m

    fx = 0.2
    fy = 0.5
    tau = 1.0

    systems = @named begin
        b1 = Body2d(m = m1, J = J1)
        b2 = Body2d(m = m2, J = J2)
        b1_b2_offset = FrameOffset(x_0 = l, y_0 = 0.0)
        b1_forces = ForceAndTorque()
        b1_fx = Constant(k = fx)
        b1_fy = Constant(k = fy)
        b1_tau = Constant(k = tau)
    end

    equations = Equation[
        connect(b1.frame, b1_forces.frame, b1_b2_offset.frame_a),
        connect(b1_forces.fx, b1_fx.output),
        connect(b1_forces.fy, b1_fy.output),
        connect(b1_forces.tau, b1_tau.output),
        connect(b2.frame, b1_b2_offset.frame_b)
    ]

    @named sys_model = System(equations, t, [], []; systems)

    sys = mtkcompile(sys_model)
    prob = ODEProblem(sys, [], (0.0, 10.0))
    sol = solve(prob)

    @test SciMLBase.successful_retcode(sol)


    # Reference implementation without the ModelingToolkitStandardLibrary for comparison
    m_total = m1 + m2
    J_total = J1 + J2 + m1 * m2 * l^2 / m_total
    l_cog = m2 / m_total * l

    ax = fx / m_total
    ay = fy / m_total

    function ode!(du, u, _, _)
        # States, with respect to the center of gravity (cog):
        # vx  ... x-velocity
        # vy  ... y-velocity
        # w   ... rotational acceleration
        # x   ... x-position
        # y   ... y-position
        # phi ... angle

        vx, vy, w, x, y, phi = u

        # Torque acting on the cog
        tau_cog = tau - l_cog * fy * cos(phi) + l_cog * fx * sin(phi)

        dvx = fx / m_total
        dvy = fy / m_total
        dw = tau_cog / J_total
        dx = vx
        dy = vy
        dphi = w

        du .= [dvx, dvy, dw, dx, dy, dphi]
    end

    u0 = [0.0, 0.0, 0.0, l_cog, 0.0, 0.0]

    sol_ref  = solve(ODEProblem(ode!, u0, (0.0, 10.0)), Tsit5())
    u_ref = sol_ref(sol.t)

    vx_cog_ref = map(x -> x[1], u_ref)
    vy_cog_ref = map(x -> x[2], u_ref)
    w_ref =      map(x -> x[3], u_ref)
    x_cog_ref =  map(x -> x[4], u_ref)
    y_cog_ref =  map(x -> x[5], u_ref)
    phi_ref =    map(x -> x[6], u_ref)

    b1_x_ref = @. x_cog_ref - l_cog * cos(phi_ref)
    b1_y_ref = @. y_cog_ref - l_cog * sin(phi_ref)
    b2_x_ref = @. x_cog_ref + (l - l_cog) * cos(phi_ref)
    b2_y_ref = @. y_cog_ref + (l - l_cog) * sin(phi_ref)

    test_tol = 1e-3
    @test all(isapprox.(sol[sys.b1.x], b1_x_ref, atol = test_tol))
    @test all(isapprox.(sol[sys.b1.y], b1_y_ref, atol = test_tol))
    @test all(isapprox.(sol[sys.b2.x], b2_x_ref, atol = test_tol))
    @test all(isapprox.(sol[sys.b2.y], b2_y_ref, atol = test_tol))
    @test all(isapprox.(sol[sys.b1.phi], phi_ref, atol = test_tol))
end

@testset "supported torque between 2d bodies" begin
    m1 = 1.0 # kg
    m2 = 0.2 # kg
    J1 = 5.0 # kg·m²
    J2 = 2.0 # kg·m²

    tau = 2.0 # N·m

    systems = @named begin
        b1 = Body2d(m = m1, J = J1, x = 0.0, y = 0.0)
        b2 = Body2d(m = m2, J = J2, x = 0.0, y = 0.0)
        supported_torque = SupportedTorque()
        b1_b2_tau = Constant(k = tau)
    end

    equations = Equation[
        connect(b1.frame, supported_torque.support),
        connect(b2.frame, supported_torque.frame),
        connect(supported_torque.tau, b1_b2_tau.output)
    ]

    @named sys_model = System(equations, t, [], []; systems)

    sys = mtkcompile(sys_model)
    prob = ODEProblem(sys, [], (0.0, 10.0))
    sol = solve(prob)

    @test SciMLBase.successful_retcode(sol)


    # Analytial solution for comparison
    f_b1_phi_ref(t) = -0.5 * tau / J1 * t^2
    f_b2_phi_ref(t) = 0.5 * tau / J2 * t^2

    b1_phi_ref = f_b1_phi_ref.(sol.t)
    b2_phi_ref = f_b2_phi_ref.(sol.t)

    test_tol = 1e-12
    @test all(isapprox.(sol[sys.b1.phi], b1_phi_ref, atol = test_tol))
    @test all(isapprox.(sol[sys.b2.phi], b2_phi_ref, atol = test_tol))
end
