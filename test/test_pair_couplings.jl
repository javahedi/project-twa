@testset "Pair couplings" begin

    @testset "Nearest-neighbor XXZ" begin
        model = SpinModel(
            Chain(5),
            XXZ(J=2.0, Δ=0.5),
        )

        @test coupling_strength(model, 1, 2) == 1.0
        @test coupling_strength(model, 1, 3) == 0.0
        @test coupling_strength(model, 2, 2) == 0.0

        @test coupling_components(model, 1, 2) == (
            2.0,
            2.0,
            1.0,
        )

        @test coupling_components(model, 1, 3) == (
            0.0,
            0.0,
            0.0,
        )

        @test coupling_components(model, 2, 2) == (
            0.0,
            0.0,
            0.0,
        )
    end


    @testset "Periodic nearest-neighbor XXZ" begin
        model = SpinModel(
            Chain(5, PeriodicBoundary()),
            XXZ(J=1.0, Δ=2.0),
        )

        @test coupling_strength(model, 1, 5) == 1.0

        @test coupling_components(model, 1, 5) == (
            1.0,
            1.0,
            2.0,
        )
    end


    @testset "Power-law XXZ" begin
        model = SpinModel(
            Chain(6),
            PowerLaw(
                XXZ(J=2.0, Δ=0.5);
                α=2.0,
            ),
        )

        @test coupling_strength(model, 1, 2) == 1.0
        @test coupling_strength(model, 1, 3) == 1 / 4
        @test coupling_strength(model, 1, 5) == 1 / 16

        @test coupling_components(model, 1, 3) == (
            0.5,
            0.5,
            0.25,
        )
    end


    @testset "Periodic power-law distance" begin
        model = SpinModel(
            Chain(6, PeriodicBoundary()),
            PowerLaw(
                XXZ(J=1.0, Δ=1.0);
                α=2.0,
            ),
        )

        # Sites 1 and 6 are nearest across the periodic boundary.
        @test coupling_strength(model, 1, 6) == 1.0

        # Sites 1 and 5 are distance 2 on the ring.
        @test coupling_strength(model, 1, 5) == 1 / 4
    end


    @testset "Symmetry" begin
        model = SpinModel(
            Chain(8, PeriodicBoundary()),
            PowerLaw(
                XXZ(J=1.3, Δ=0.7);
                α=3.0,
            ),
        )

        @test coupling_strength(model, 2, 7) ==
              coupling_strength(model, 7, 2)

        @test coupling_components(model, 2, 7) ==
              coupling_components(model, 7, 2)
    end


    @testset "Numeric types" begin
        model = SpinModel(
            Chain(4),
            XXZ(J=Float32(1), Δ=Float32(0.5)),
        )

        @test coupling_strength(model, 1, 2) isa Float32

        components = coupling_components(model, 1, 2)

        @test components isa NTuple{3,Float32}
        @test components == (1.0f0, 1.0f0, 0.5f0)
    end


    @testset "Bounds checking" begin
        model = SpinModel(
            Chain(5),
            XXZ(),
        )

        @test_throws BoundsError coupling_strength(model, 0, 1)
        @test_throws BoundsError coupling_strength(model, 1, 6)

        @test_throws BoundsError coupling_components(model, 0, 1)
        @test_throws BoundsError coupling_components(model, 1, 6)
    end
end



@testset "Ising interactions" begin
    @testset "Nearest-neighbor axes" begin
        mx = SpinModel(
            Chain(4),
            Ising(:x; J=2.0),
        )

        my = SpinModel(
            Chain(4),
            Ising(:y; J=2.0),
        )

        mz = SpinModel(
            Chain(4),
            Ising(:z; J=2.0),
        )

        @test coupling_components(mx, 1, 2) ==
            (2.0, 0.0, 0.0)

        @test coupling_components(my, 1, 2) ==
            (0.0, 2.0, 0.0)

        @test coupling_components(mz, 1, 2) ==
            (0.0, 0.0, 2.0)

        @test coupling_components(mx, 1, 3) ==
            (0.0, 0.0, 0.0)
    end

    @testset "Power-law Ising" begin
        model = SpinModel(
            Chain(5),
            PowerLaw(
                Ising(:x; J=2.0);
                α=2.0,
            ),
        )

        components_13 = coupling_components(
            model,
            1,
            3,
        )

        @test components_13[1] ≈ 0.5
        @test components_13[2] == 0.0
        @test components_13[3] == 0.0
    end

    @testset "Validation" begin
        @test_throws ArgumentError Ising(:q; J=1.0)
    end
end