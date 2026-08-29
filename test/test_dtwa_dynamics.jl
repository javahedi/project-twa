@testset "DTWA dynamics" begin

    @testset "Zero field and aligned isotropic state" begin
        model = SpinModel(
            Chain(3),
            XXZ(J=1.0, Δ=1.0),
        )

        u = [
            0.0 0.0 1.0
            0.0 0.0 1.0
            0.0 0.0 1.0
        ]

        du = similar(u)

        dtwa_rhs!(
            du,
            u,
            model,
            0.0,
        )

        @test du == zeros(3, 3)
    end


    @testset "Single-spin field precession" begin
        model = SpinModel(
            Chain(1),
            Field(:z, 0.5),
        )

        u = reshape(
            [1.0, 0.0, 0.0],
            1,
            3,
        )

        du = similar(u)

        dtwa_rhs!(
            du,
            u,
            model,
            0.0,
        )

        # H = 0.5 σᶻ gives ds/dt = 2 h × s,
        # hence +ŷ for a spin initially along +x.
        @test du ≈ reshape(
            [0.0, 1.0, 0.0],
            1,
            3,
        )
    end


    @testset "Multiple fields add" begin
        model = SpinModel(
            Chain(1),
            Field(:x, 0.25) +
            Field(:z, 0.5),
        )

        u = reshape(
            [0.0, 1.0, 0.0],
            1,
            3,
        )

        du = similar(u)

        dtwa_rhs!(
            du,
            u,
            model,
            0.0,
        )

        @test du ≈ reshape(
            [-1.0, 0.0, 0.5],
            1,
            3,
        )
    end


    @testset "Two-spin XXZ interaction" begin
        model = SpinModel(
            Chain(2),
            XXZ(J=1.0, Δ=0.5),
        )

        u = [
            1.0 0.0 0.0
            0.0 0.0 1.0
        ]

        du = similar(u)

        dtwa_rhs!(
            du,
            u,
            model,
            0.0,
        )

        # Site 1 sees h_eff = (0, 0, 0.5).
        # Site 2 sees h_eff = (1, 0, 0).
        @test du ≈ [
             0.0  1.0 0.0
             0.0 -2.0 0.0
        ]
    end


    @testset "Power-law interaction participates" begin
        model = SpinModel(
            Chain(3),
            PowerLaw(
                XXZ(J=2.0, Δ=1.0);
                α=1.0,
            ),
        )

        u = [
            1.0 0.0 0.0
            0.0 0.0 0.0
            0.0 0.0 1.0
        ]

        du = similar(u)

        dtwa_rhs!(
            du,
            u,
            model,
            0.0,
        )

        # Distance between sites 1 and 3 is 2, so J/r = 1.
        @test du[1, :] ≈ [0.0, 2.0, 0.0]
        @test du[3, :] ≈ [0.0, -2.0, 0.0]
    end


    @testset "Instantaneous spin length is conserved" begin
        model = SpinModel(
            Chain(
                4,
                PeriodicBoundary(),
            ),
            XXZ(J=0.7, Δ=1.3) +
            Field(:x, 0.2) +
            Field(:z, -0.4),
        )

        u = [
             1.0  0.2 -0.4
            -0.3  1.1  0.5
             0.7 -0.8  1.2
            -1.0  0.4  0.6
        ]

        du = similar(u)

        dtwa_rhs!(
            du,
            u,
            model,
            0.0,
        )

        # d|sᵢ|²/dt = 2 sᵢ⋅dsᵢ/dt = 0 for precession.
        for i in axes(u, 1)
            tangent = (
                u[i, 1] * du[i, 1] +
                u[i, 2] * du[i, 2] +
                u[i, 3] * du[i, 3]
            )

            @test isapprox(
                tangent,
                0.0;
                atol=1e-12,
            )
        end
    end


    @testset "Float32 trajectory" begin
        model = SpinModel(
            Chain(2),
            XXZ(
                J=Float32(1),
                Δ=Float32(0.5),
            ),
        )

        u = Float32[
            1 0 0
            0 0 1
        ]

        du = similar(u)

        dtwa_rhs!(
            du,
            u,
            model,
            0.0f0,
        )

        @test eltype(du) === Float32
        @test du == Float32[
             0  1 0
             0 -2 0
        ]
    end


    @testset "Layout validation" begin
        model = SpinModel(
            Chain(3),
            XXZ(),
        )

        u_bad = zeros(3, 2)
        du_bad = similar(u_bad)

        @test_throws DimensionMismatch dtwa_rhs!(
            du_bad,
            u_bad,
            model,
            0.0,
        )

        u = zeros(3, 3)
        du = zeros(2, 3)

        @test_throws DimensionMismatch dtwa_rhs!(
            du,
            u,
            model,
            0.0,
        )
    end
end
