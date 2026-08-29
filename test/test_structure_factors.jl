@testset "Structure factors" begin

    @testset "Polarized state on default momentum grid" begin
        model = SpinModel(
            Chain(
                4,
                PeriodicBoundary(),
            ),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=8);
            tspan=(0.0, 1.0),
            saveat=[0.0, 1.0],
            rng=Xoshiro(40),
        )

        sf = structure_factor(
            result,
            geometry(model),
            :z,
        )

        @test sf isa StructureFactor
        @test length(sf.momenta) == 4
        @test size(sf.values) == (2, 4)
        @test sf.axis == :z
        @test sf.connected == false

        @test sf.momenta ≈ [
            0.0,
            pi / 2,
            pi,
            3pi / 2,
        ]

        @test sf.values[:, 1] ≈ [4.0, 4.0]
        @test sf.values[:, 2] ≈ [0.0, 0.0] atol=1e-12
        @test sf.values[:, 3] ≈ [0.0, 0.0] atol=1e-12
        @test sf.values[:, 4] ≈ [0.0, 0.0] atol=1e-12
    end


    @testset "Connected polarized structure factor vanishes" begin
        model = SpinModel(
            Chain(5),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=8);
            tspan=(0.0, 0.5),
            saveat=[0.0, 0.5],
            rng=Xoshiro(41),
        )

        sf = structure_factor(
            result,
            geometry(model),
            :z;
            connected=true,
        )

        @test sf.connected == true
        @test all(
            isapprox.(sf.values, 0.0; atol=1e-12)
        )
    end


    @testset "Custom momenta" begin
        model = SpinModel(
            Chain(3),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=4);
            tspan=(0.0, 0.2),
            saveat=[0.0, 0.2],
            rng=Xoshiro(42),
        )

        sf = structure_factor(
            result,
            geometry(model),
            :z;
            momenta=[0.0, pi],
        )

        @test sf.momenta ≈ [0.0, pi]
        @test size(sf.values) == (2, 2)

        @test sf.values[:, 1] ≈ [3.0, 3.0]

        # For three identical spins:
        #
        # (1/N) |1 - 1 + 1|² = 1/3.
        @test sf.values[:, 2] ≈ [1 / 3, 1 / 3]
    end


    @testset "Independent transverse fluctuations" begin
        # Exact balanced ensemble for two z-polarized spins:
        #
        # <σx_i> = 0,
        # <σx_1 σx_2> = 0,
        # <(σx_i)^2> = 1.
        #
        # Therefore both ordinary and connected Sx(q) equal 1 for all q.
        trajectories = zeros(Float64, 2, 3, 1, 4)

        x_samples = (
            (1.0, 1.0),
            (1.0, -1.0),
            (-1.0, 1.0),
            (-1.0, -1.0),
        )

        for k in 1:4
            trajectories[1, 1, 1, k] = x_samples[k][1]
            trajectories[2, 1, 1, k] = x_samples[k][2]

            trajectories[1, 3, 1, k] = 1.0
            trajectories[2, 3, 1, k] = 1.0
        end

        result = DTWAResult(
            [0.0],
            trajectories,
        )

        sf = structure_factor(
            result,
            Chain(2),
            :x,
        )

        sf_connected = structure_factor(
            result,
            Chain(2),
            :x;
            connected=true,
        )

        @test sf.values ≈ ones(1, 2)
        @test sf_connected.values ≈ ones(1, 2)
    end


    @testset "q=0 relation to collective second moment" begin
        model = SpinModel(
            Chain(4),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            DomainWall(),
            DTWA(trajectories=8);
            tspan=(0.0, 1.0),
            saveat=[0.0, 1.0],
            rng=Xoshiro(43),
        )

        sf = structure_factor(
            result,
            geometry(model),
            :z;
            momenta=[0.0],
        )

        second = expectation(
            result,
            SpinSecondMoment(:z),
        )

        # With our conventions:
        #
        # S(q=0) = <S²>/N
        #        = N * [<S²>/N²].
        @test sf.values[:, 1] ≈
            geometry(model).nsites .* second
    end


    @testset "Connected q=0 relation to collective variance" begin
        trajectories = zeros(Float64, 2, 3, 1, 4)

        x_samples = (
            (1.0, 1.0),
            (1.0, -1.0),
            (-1.0, 1.0),
            (-1.0, -1.0),
        )

        for k in 1:4
            trajectories[1, 1, 1, k] = x_samples[k][1]
            trajectories[2, 1, 1, k] = x_samples[k][2]
        end

        result = DTWAResult(
            [0.0],
            trajectories,
        )

        sf = structure_factor(
            result,
            Chain(2),
            :x;
            connected=true,
            momenta=[0.0],
        )

        variance = expectation(
            result,
            SpinVariance(:x),
        )

        @test sf.values[:, 1] ≈
            2 .* variance
    end


    @testset "Geometry/result mismatch" begin
        trajectories = zeros(Float64, 2, 3, 1, 1)

        result = DTWAResult(
            [0.0],
            trajectories,
        )

        @test_throws DimensionMismatch structure_factor(
            result,
            Chain(3),
            :z,
        )
    end


    @testset "Input validation" begin
        trajectories = zeros(Float64, 1, 3, 1, 1)

        result = DTWAResult(
            [0.0],
            trajectories,
        )

        @test_throws ArgumentError structure_factor(
            result,
            Chain(1),
            :q,
        )

        @test_throws ArgumentError structure_factor(
            result,
            Chain(1),
            :z;
            momenta=Float64[],
        )
    end
end
