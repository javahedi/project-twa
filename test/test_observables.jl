using Random

@testset "Observables" begin

    @testset "Observable construction" begin
        @test LocalMagnetization(:x, 2) isa AbstractObservable
        @test Magnetization(:z) isa AbstractObservable
        @test Correlation(:x, 1, :z, 2) isa AbstractObservable
        @test ConnectedCorrelation(:x, 1, :x, 2) isa AbstractObservable

        @test_throws ArgumentError LocalMagnetization(:q, 1)
        @test_throws ArgumentError Magnetization(:q)
        @test_throws ArgumentError LocalMagnetization(:x, 0)
        @test_throws ArgumentError Correlation(:x, 1, :y, 1)
    end


    @testset "Polarized-state magnetization" begin
        model = SpinModel(
            Chain(4),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=64);
            tspan=(0.0, 1.0),
            saveat=[0.0, 1.0],
            rng=Xoshiro(1),
        )

        @test expectation(
            result,
            LocalMagnetization(:z, 3),
        ) == [1.0, 1.0]

        @test expectation(
            result,
            Magnetization(:z),
        ) == [1.0, 1.0]
    end


    @testset "Domain-wall magnetization density" begin
        model = SpinModel(
            Chain(4),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            DomainWall(),
            DTWA(trajectories=32);
            tspan=(0.0, 1.0),
            saveat=[0.0, 1.0],
            rng=Xoshiro(2),
        )

        @test expectation(
            result,
            Magnetization(:z),
        ) == [0.0, 0.0]

        @test expectation(
            result,
            LocalMagnetization(:z, 1),
        ) == [1.0, 1.0]

        @test expectation(
            result,
            LocalMagnetization(:z, 4),
        ) == [-1.0, -1.0]
    end


    @testset "Distinct-site correlations" begin
        model = SpinModel(
            Chain(3),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=128);
            tspan=(0.0, 1.0),
            saveat=[0.0, 1.0],
            rng=Xoshiro(3),
        )

        @test expectation(
            result,
            Correlation(:z, 1, :z, 3),
        ) == [1.0, 1.0]

        @test expectation(
            result,
            ConnectedCorrelation(:z, 1, :z, 3),
        ) == [0.0, 0.0]
    end


    @testset "Same-site Pauli identity" begin
        model = SpinModel(
            Chain(2),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Polarized(:x),
            DTWA(trajectories=16);
            tspan=(0.0, 0.5),
            saveat=[0.0, 0.5],
            rng=Xoshiro(4),
        )

        @test expectation(
            result,
            Correlation(:x, 1, :x, 1),
        ) == [1.0, 1.0]

        @test expectation(
            result,
            Correlation(:y, 2, :y, 2),
        ) == [1.0, 1.0]
    end


    @testset "Connected same-site variance identity" begin
        model = SpinModel(
            Chain(1),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=32);
            tspan=(0.0, 1.0),
            saveat=[0.0, 1.0],
            rng=Xoshiro(5),
        )

        @test expectation(
            result,
            ConnectedCorrelation(:z, 1, :z, 1),
        ) == [0.0, 0.0]
    end


    @testset "Observable output follows saved times" begin
        model = SpinModel(
            Chain(2),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Down(),
            DTWA(trajectories=8);
            tspan=(0.0, 1.0),
            saveat=0.2,
            rng=Xoshiro(6),
        )

        values = expectation(
            result,
            Magnetization(:z),
        )

        @test length(values) == length(times(result))
        @test all(values .== -1.0)
    end


    @testset "Site validation" begin
        model = SpinModel(
            Chain(2),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=4);
            tspan=(0.0, 1.0),
            saveat=[0.0, 1.0],
            rng=Xoshiro(7),
        )

        @test_throws BoundsError expectation(
            result,
            LocalMagnetization(:z, 3),
        )

        @test_throws BoundsError expectation(
            result,
            Correlation(:z, 1, :z, 3),
        )
    end
end
