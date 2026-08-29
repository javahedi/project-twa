@testset "Compiled DTWA dynamics" begin

    @testset "Compilation stores model size and field" begin
        model = SpinModel(
            Chain(4),
            XXZ(J=1.0, Δ=0.5) +
            Field(:x, 0.2) +
            Field(:z, -0.3),
        )

        compiled = compile(
            model,
            DTWA(trajectories=10),
        )

        @test compiled isa CompiledDTWA
        @test compiled.nsites == 4
        @test compiled.field == (0.2, 0.0, -0.3)
    end


    @testset "Nearest-neighbor compilation stores only bonds" begin
        model = SpinModel(
            Chain(6),
            XXZ(J=1.0, Δ=0.5),
        )

        compiled = compile(
            model,
            DTWA(),
        )

        @test length(compiled.pairs) == 5

        @test compiled.pairs[1] ==
            (1, 2, 1.0, 1.0, 0.5)

        @test compiled.pairs[end] ==
            (5, 6, 1.0, 1.0, 0.5)
    end


    @testset "Periodic nearest-neighbor compilation" begin
        model = SpinModel(
            Chain(
                5,
                PeriodicBoundary(),
            ),
            XXZ(J=2.0, Δ=0.25),
        )

        compiled = compile(
            model,
            DTWA(),
        )

        @test length(compiled.pairs) == 5

        pair_sites = Set(
            (pair[1], pair[2])
            for pair in compiled.pairs
        )

        @test pair_sites == Set([
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (1, 5),
        ])
    end


    @testset "Power law stores all nonzero pairs" begin
        model = SpinModel(
            Chain(4),
            PowerLaw(
                XXZ(J=1.0, Δ=1.0);
                α=2.0,
            ),
        )

        compiled = compile(
            model,
            DTWA(),
        )

        @test length(compiled.pairs) == 6

        @test compiled.pairs[1][1:2] == (1, 2)

        @test compiled.pairs[end][1:2] == (3, 4)
    end


    @testset "Reference and compiled RHS agree" begin
        models = (
            SpinModel(
                Chain(5),
                XXZ(J=0.7, Δ=1.2),
            ),
            SpinModel(
                Chain(
                    5,
                    PeriodicBoundary(),
                ),
                XXZ(J=0.4, Δ=0.6) +
                Field(:y, 0.3),
            ),
            SpinModel(
                Chain(5),
                PowerLaw(
                    XXZ(J=0.9, Δ=1.1);
                    α=1.5,
                ) +
                Field(:x, -0.2) +
                Field(:z, 0.1),
            ),
        )

        u = [
             1.0  0.2 -0.4
            -0.3  1.1  0.5
             0.7 -0.8  1.2
            -1.0  0.4  0.6
             0.1 -0.5  0.9
        ]

        for model in models
            compiled = compile(
                model,
                DTWA(),
            )

            du_reference = similar(u)
            du_compiled = similar(u)

            dtwa_rhs!(
                du_reference,
                u,
                model,
                0.0,
            )

            dtwa_rhs!(
                du_compiled,
                u,
                compiled,
                0.0,
            )

            @test du_compiled ≈ du_reference
        end
    end


    @testset "Float32 compilation" begin
        model = SpinModel(
            Chain(3),
            XXZ(
                J=Float32(1),
                Δ=Float32(0.5),
            ) +
            Field(
                :z,
                Float32(0.25),
            ),
        )

        compiled = compile(
            model,
            DTWA(),
        )

        @test compiled.field isa NTuple{3,Float32}

        @test all(
            pair ->
                pair[3] isa Float32 &&
                pair[4] isa Float32 &&
                pair[5] isa Float32,
            compiled.pairs,
        )

        u = Float32[
            1 0 0
            0 1 0
            0 0 1
        ]

        du = similar(u)

        dtwa_rhs!(
            du,
            u,
            compiled,
            0.0f0,
        )

        @test eltype(du) === Float32
    end


    @testset "Compiled spin-length conservation" begin
        model = SpinModel(
            Chain(
                4,
                PeriodicBoundary(),
            ),
            XXZ(J=0.6, Δ=1.4) +
            Field(:x, 0.25) +
            Field(:z, -0.35),
        )

        compiled = compile(
            model,
            DTWA(),
        )

        u = [
             0.9  0.2 -0.3
            -0.4  1.0  0.6
             0.8 -0.7  1.1
            -0.9  0.5  0.4
        ]

        du = similar(u)

        dtwa_rhs!(
            du,
            u,
            compiled,
            0.0,
        )

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


    @testset "Compiled layout validation" begin
        model = SpinModel(
            Chain(3),
            XXZ(),
        )

        compiled = compile(
            model,
            DTWA(),
        )

        u = zeros(2, 3)
        du = similar(u)

        @test_throws DimensionMismatch dtwa_rhs!(
            du,
            u,
            compiled,
            0.0,
        )
    end
end
