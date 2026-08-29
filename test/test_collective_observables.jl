@testset "Collective observables" begin

    @testset "Construction" begin
        @test SpinSecondMoment(:x) isa AbstractObservable
        @test SpinVariance(:z) isa AbstractObservable

        @test_throws ArgumentError SpinSecondMoment(:q)
        @test_throws ArgumentError SpinVariance(:q)
    end


    @testset "Fully polarized collective spin" begin
        model = SpinModel(
            Chain(4),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=16);
            tspan=(0.0, 1.0),
            saveat=[0.0, 1.0],
            rng=Xoshiro(20),
        )

        @test expectation(
            result,
            SpinSecondMoment(:z),
        ) == [1.0, 1.0]

        @test expectation(
            result,
            SpinVariance(:z),
        ) == [0.0, 0.0]
    end


    @testset "Domain wall has zero longitudinal collective spin" begin
        model = SpinModel(
            Chain(4),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            DomainWall(),
            DTWA(trajectories=16);
            tspan=(0.0, 1.0),
            saveat=[0.0, 1.0],
            rng=Xoshiro(21),
        )

        @test expectation(
            result,
            SpinSecondMoment(:z),
        ) == [0.0, 0.0]

        @test expectation(
            result,
            SpinVariance(:z),
        ) == [0.0, 0.0]
    end


    @testset "Transverse product-state variance" begin
        # Exact four-trajectory representation of the independent DTWA
        # transverse samples for two spins polarized along z:
        #
        #     (x₁, x₂) = (++), (+-), (-+), (--).
        #
        # Therefore
        #
        #     ⟨Sx⟩ = 0,
        #     ⟨Sx²⟩ = 2,
        #
        # and, for N = 2,
        #
        #     ⟨Sx²⟩ / N² = (ΔSx)² / N² = 1/2.
        #
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

        @test expectation(
            result,
            Magnetization(:x),
        ) == [0.0]

        @test expectation(
            result,
            Correlation(:x, 1, :x, 2),
        ) == [0.0]

        @test expectation(
            result,
            SpinSecondMoment(:x),
        ) == [0.5]

        @test expectation(
            result,
            SpinVariance(:x),
        ) == [0.5]
    end


    @testset "Variance identity" begin
        model = SpinModel(
            Chain(3),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=32);
            tspan=(0.0, 0.5),
            saveat=[0.0, 0.5],
            rng=Xoshiro(22),
        )

        second = expectation(
            result,
            SpinSecondMoment(:z),
        )

        mean = expectation(
            result,
            Magnetization(:z),
        )

        variance = expectation(
            result,
            SpinVariance(:z),
        )

        @test variance ≈ second .- mean .^ 2
    end


    @testset "Single-site limit" begin
        model = SpinModel(
            Chain(1),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Polarized(:x),
            DTWA(trajectories=8);
            tspan=(0.0, 1.0),
            saveat=[0.0, 1.0],
            rng=Xoshiro(23),
        )

        @test expectation(
            result,
            SpinSecondMoment(:x),
        ) == [1.0, 1.0]

        @test expectation(
            result,
            SpinVariance(:x),
        ) == [0.0, 0.0]
    end
end
