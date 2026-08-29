using Random

@testset "DTWA initial ensemble" begin

    @testset "Shape" begin
        model = SpinModel(
            Chain(6),
            XXZ(),
        )

        method = DTWA(
            trajectories=12,
        )

        samples = sample_initial(
            model,
            Up(),
            method;
            rng=Xoshiro(1),
        )

        @test size(samples) == (6, 3, 12)
        @test eltype(samples) === Float64
    end


    @testset "Float32 support" begin
        model = SpinModel(
            Chain(4),
            XXZ(),
        )

        samples = sample_initial(
            model,
            Down(),
            DTWA(trajectories=5);
            rng=Xoshiro(2),
            T=Float32,
        )

        @test size(samples) == (4, 3, 5)
        @test eltype(samples) === Float32
    end


    @testset "Up-state longitudinal component" begin
        model = SpinModel(
            Chain(8),
            XXZ(),
        )

        samples = sample_initial(
            model,
            Up(),
            DTWA(trajectories=20);
            rng=Xoshiro(3),
        )

        @test all(samples[:, 3, :] .== 1)

        @test all(
            value -> value in (-1, 1),
            samples[:, 1, :],
        )

        @test all(
            value -> value in (-1, 1),
            samples[:, 2, :],
        )
    end


    @testset "Domain-wall structure" begin
        model = SpinModel(
            Chain(6),
            XXZ(),
        )

        samples = sample_initial(
            model,
            DomainWall(),
            DTWA(trajectories=10);
            rng=Xoshiro(4),
        )

        expected = [
            1,
            1,
            1,
            -1,
            -1,
            -1,
        ]

        for trajectory in axes(samples, 3)
            @test samples[:, 3, trajectory] == expected
        end
    end


    @testset "Trajectory slices preserve N×3 layout" begin
        model = SpinModel(
            Chain(7),
            XXZ(),
        )

        samples = sample_initial(
            model,
            Polarized(:x),
            DTWA(trajectories=4);
            rng=Xoshiro(5),
        )

        trajectory = @view samples[:, :, 2]

        @test size(trajectory) == (7, 3)
        @test all(trajectory[:, 1] .== 1)
    end


    @testset "Explicit RNG is reproducible" begin
        model = SpinModel(
            Chain(10),
            XXZ(),
        )

        method = DTWA(
            trajectories=30,
        )

        samples_a = sample_initial(
            model,
            Up(),
            method;
            rng=Xoshiro(1234),
        )

        samples_b = sample_initial(
            model,
            Up(),
            method;
            rng=Xoshiro(1234),
        )

        @test samples_a == samples_b
    end


    @testset "Different trajectories are sampled independently" begin
        model = SpinModel(
            Chain(20),
            XXZ(),
        )

        samples = sample_initial(
            model,
            Up(),
            DTWA(trajectories=20);
            rng=Xoshiro(99),
        )

        first = @view samples[:, :, 1]

        @test any(
            trajectory -> samples[:, :, trajectory] != first,
            2:size(samples, 3),
        )
    end


    @testset "Periodic geometry uses the same site count" begin
        model = SpinModel(
            Chain(
                5,
                PeriodicBoundary(),
            ),
            XXZ(),
        )

        samples = sample_initial(
            model,
            Down(),
            DTWA(trajectories=3);
            rng=Xoshiro(7),
        )

        @test size(samples) == (5, 3, 3)
    end
end
