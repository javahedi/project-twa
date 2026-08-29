@testset "Correlation profiles" begin

    @testset "Open-chain distances" begin
        model = SpinModel(
            Chain(4),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=8);
            tspan=(0.0, 1.0),
            saveat=[0.0, 1.0],
            rng=Xoshiro(30),
        )

        profile = correlation_profile(
            result,
            geometry(model),
            :z,
        )

        @test profile isa CorrelationProfile
        @test profile.distances == [0, 1, 2, 3]
        @test size(profile.values) == (2, 4)
        @test profile.axis == :z
        @test profile.connected == false
        @test all(profile.values .== 1.0)
    end


    @testset "Periodic-chain distances" begin
        model = SpinModel(
            Chain(
                6,
                PeriodicBoundary(),
            ),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=8);
            tspan=(0.0, 0.5),
            saveat=[0.0, 0.5],
            rng=Xoshiro(31),
        )

        profile = correlation_profile(
            result,
            geometry(model),
            :z,
        )

        @test profile.distances == [0, 1, 2, 3]
        @test all(profile.values .== 1.0)
    end


    @testset "Odd periodic maximum distance" begin
        model = SpinModel(
            Chain(
                5,
                PeriodicBoundary(),
            ),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=4);
            tspan=(0.0, 0.2),
            saveat=[0.0, 0.2],
            rng=Xoshiro(32),
        )

        profile = correlation_profile(
            result,
            geometry(model),
            :z,
        )

        @test profile.distances == [0, 1, 2]
    end


    @testset "Domain-wall spatial profile" begin
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
            rng=Xoshiro(33),
        )

        profile = correlation_profile(
            result,
            geometry(model),
            :z,
        )

        # For (+,+,-,-):
        #
        # r = 0:  1
        # r = 1: (1 - 1 + 1) / 3 = 1/3
        # r = 2: (-1 - 1) / 2 = -1
        # r = 3: -1
        expected = [
            1.0,
            1 / 3,
            -1.0,
            -1.0,
        ]

        @test profile.values[1, :] ≈ expected
        @test profile.values[2, :] ≈ expected
    end


    @testset "Connected profile for deterministic longitudinal state" begin
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
            rng=Xoshiro(34),
        )

        profile = correlation_profile(
            result,
            geometry(model),
            :z;
            connected=true,
        )

        @test profile.connected == true
        @test all(profile.values .== 0.0)
    end


    @testset "On-site transverse connected correlation" begin
        # Construct an exactly balanced one-site ensemble with <σx> = 0.
        trajectories = zeros(Float64, 1, 3, 1, 2)

        trajectories[1, 1, 1, 1] = 1.0
        trajectories[1, 1, 1, 2] = -1.0

        result = DTWAResult(
            [0.0],
            trajectories,
        )

        profile = correlation_profile(
            result,
            Chain(1),
            :x;
            connected=true,
        )

        @test profile.distances == [0]
        @test profile.values == reshape([1.0], 1, 1)
    end


    @testset "Periodic pair counting" begin
        chain = Chain(
            6,
            PeriodicBoundary(),
        )

        @test length(
            TWA._site_pairs_at_distance(
                chain,
                1,
            ),
        ) == 6

        @test length(
            TWA._site_pairs_at_distance(
                chain,
                2,
            ),
        ) == 6

        # Opposite points on an even ring occur only once per unordered pair.
        @test length(
            TWA._site_pairs_at_distance(
                chain,
                3,
            ),
        ) == 3
    end


    @testset "Geometry/result mismatch" begin
        model = SpinModel(
            Chain(3),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=4);
            tspan=(0.0, 1.0),
            saveat=[0.0, 1.0],
            rng=Xoshiro(35),
        )

        @test_throws DimensionMismatch correlation_profile(
            result,
            Chain(4),
            :z,
        )
    end


    @testset "Axis validation" begin
        trajectories = zeros(Float64, 1, 3, 1, 1)

        result = DTWAResult(
            [0.0],
            trajectories,
        )

        @test_throws ArgumentError correlation_profile(
            result,
            Chain(1),
            :q,
        )
    end
end
