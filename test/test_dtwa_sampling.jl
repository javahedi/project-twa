using Random

@testset "DTWA sampling" begin

    @testset "Shape and element type" begin
        sample = sample_dtwa(
            Up(),
            6;
            rng=Xoshiro(1),
        )

        @test size(sample) == (6, 3)
        @test eltype(sample) === Float64

        sample32 = sample_dtwa(
            Up(),
            6;
            rng=Xoshiro(1),
            T=Float32,
        )

        @test eltype(sample32) === Float32
    end


    @testset "Up state" begin
        sample = sample_dtwa(
            Up(),
            20;
            rng=Xoshiro(2),
        )

        @test all(sample[:, 3] .== 1)

        @test all(
            value -> value in (-1, 1),
            sample[:, 1],
        )

        @test all(
            value -> value in (-1, 1),
            sample[:, 2],
        )
    end


    @testset "Down state" begin
        sample = sample_dtwa(
            Down(),
            20;
            rng=Xoshiro(3),
        )

        @test all(sample[:, 3] .== -1)

        @test all(
            value -> value in (-1, 1),
            sample[:, 1],
        )

        @test all(
            value -> value in (-1, 1),
            sample[:, 2],
        )
    end


    @testset "Polarized x state" begin
        sample = sample_dtwa(
            Polarized(:x; sign=-1),
            20;
            rng=Xoshiro(4),
        )

        @test all(sample[:, 1] .== -1)

        @test all(
            value -> value in (-1, 1),
            sample[:, 2],
        )

        @test all(
            value -> value in (-1, 1),
            sample[:, 3],
        )
    end


    @testset "Domain wall" begin
        sample = sample_dtwa(
            DomainWall(),
            6;
            rng=Xoshiro(5),
        )

        @test sample[:, 3] == [
            1,
            1,
            1,
            -1,
            -1,
            -1,
        ]

        @test all(
            value -> value in (-1, 1),
            sample[:, 1],
        )

        @test all(
            value -> value in (-1, 1),
            sample[:, 2],
        )
    end


    @testset "Odd domain wall" begin
        sample = sample_dtwa(
            DomainWall(axis=:y),
            5;
            rng=Xoshiro(6),
        )

        @test sample[:, 2] == [
            1,
            1,
            1,
            -1,
            -1,
        ]
    end


    @testset "Chain convenience overload" begin
        chain = Chain(
            8,
            PeriodicBoundary(),
        )

        sample = sample_dtwa(
            Up(),
            chain;
            rng=Xoshiro(7),
        )

        @test size(sample) == (8, 3)
        @test all(sample[:, 3] .== 1)
    end


    @testset "Reproducibility with explicit RNG" begin
        sample_a = sample_dtwa(
            Up(),
            30;
            rng=Xoshiro(1234),
        )

        sample_b = sample_dtwa(
            Up(),
            30;
            rng=Xoshiro(1234),
        )

        @test sample_a == sample_b
    end


    @testset "Independent transverse sampling" begin
        sample = sample_dtwa(
            Up(),
            200;
            rng=Xoshiro(99),
        )

        transverse_pairs = Set(
            (sample[i, 1], sample[i, 2])
            for i in axes(sample, 1)
        )

        @test transverse_pairs == Set([
            (-1.0, -1.0),
            (-1.0,  1.0),
            ( 1.0, -1.0),
            ( 1.0,  1.0),
        ])
    end


    @testset "Validation" begin
        @test_throws ArgumentError sample_dtwa(
            Up(),
            0,
        )

        @test_throws ArgumentError sample_dtwa(
            Up(),
            -3,
        )
    end
end
