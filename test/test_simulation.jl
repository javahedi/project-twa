using Random

@testset "Simulation API" begin

    @testset "Basic DTWA simulation" begin
        model = SpinModel(
            Chain(3),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=4);
            tspan=(0.0, 1.0),
            saveat=0.25,
            rng=Xoshiro(1),
        )

        @test result isa DTWAResult
        @test times(result) == [
            0.0,
            0.25,
            0.5,
            0.75,
            1.0,
        ]

        @test size(result.trajectories) ==
            (3, 3, 5, 4)

        @test ntrajectories(result) == 4
    end


    @testset "Zero Hamiltonian dynamics are static" begin
        model = SpinModel(
            Chain(5),
            Field(:x, 0.0),
        )

        result = simulate(
            model,
            DomainWall(),
            DTWA(trajectories=3);
            tspan=(0.0, 2.0),
            saveat=[0.0, 0.5, 1.0, 2.0],
            rng=Xoshiro(12),
        )

        for k in 1:ntrajectories(result)
            first_state =
                @view result.trajectories[:, :, 1, k]

            for time_index in axes(
                result.trajectories,
                3,
            )
                @test result.trajectories[
                    :,
                    :,
                    time_index,
                    k,
                ] == first_state
            end
        end
    end


    @testset "Single-spin precession preserves z" begin
        model = SpinModel(
            Chain(1),
            Field(:z, 0.5),
        )

        result = simulate(
            model,
            Polarized(:x),
            DTWA(trajectories=5);
            tspan=(0.0, 1.0),
            saveat=range(
                0.0,
                1.0;
                length=6,
            ),
            rng=Xoshiro(5),
            abstol=1e-10,
            reltol=1e-10,
        )

        # Under a z-directed field, z is conserved for every
        # independently sampled trajectory.
        for k in 1:ntrajectories(result)
            z0 = result.trajectories[
                1,
                3,
                1,
                k,
            ]

            @test all(
                result.trajectories[
                    1,
                    3,
                    :,
                    k,
                ] .≈ z0
            )
        end
    end


    @testset "Spin length is preserved through integration" begin
        model = SpinModel(
            Chain(1),
            Field(:z, 0.7),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=2);
            tspan=(0.0, 3.0),
            saveat=0.2,
            rng=Xoshiro(9),
            abstol=1e-11,
            reltol=1e-11,
        )

        for k in 1:ntrajectories(result)
            initial_norm2 = sum(
                abs2,
                @view result.trajectories[
                    1,
                    :,
                    1,
                    k,
                ]
            )

            for time_index in axes(
                result.trajectories,
                3,
            )
                norm2 = sum(
                    abs2,
                    @view result.trajectories[
                        1,
                        :,
                        time_index,
                        k,
                    ]
                )

                @test isapprox(
                    norm2,
                    initial_norm2;
                    atol=1e-8,
                )
            end
        end
    end


    @testset "Default output grid" begin
        model = SpinModel(
            Chain(1),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=1);
            tspan=(0.0, 4.0),
            rng=Xoshiro(3),
        )

        @test length(times(result)) == 101
        @test first(times(result)) == 0.0
        @test last(times(result)) == 4.0
    end


    @testset "Float32 simulation" begin
        model = SpinModel(
            Chain(1),
            Field(
                :z,
                Float32(0.5),
            ),
        )

        result = simulate(
            model,
            Up(),
            DTWA(trajectories=2);
            tspan=(0.0f0, 1.0f0),
            saveat=0.25f0,
            rng=Xoshiro(4),
            T=Float32,
        )

        @test eltype(result.t) === Float32
        @test eltype(result.trajectories) === Float32
    end


    @testset "Trajectory accessor" begin
        model = SpinModel(
            Chain(2),
            Field(:z, 0.0),
        )

        result = simulate(
            model,
            Down(),
            DTWA(trajectories=3);
            tspan=(0.0, 1.0),
            saveat=[0.0, 1.0],
            rng=Xoshiro(8),
        )

        tr = trajectory(result, 2)

        @test size(tr) == (2, 3, 2)

        @test_throws BoundsError trajectory(
            result,
            0,
        )

        @test_throws BoundsError trajectory(
            result,
            4,
        )
    end


    @testset "Input validation" begin
        model = SpinModel(
            Chain(1),
            Field(:z, 0.0),
        )

        method = DTWA(
            trajectories=1,
        )

        @test_throws ArgumentError simulate(
            model,
            Up(),
            method;
            tspan=(1.0, 0.0),
        )

        @test_throws ArgumentError simulate(
            model,
            Up(),
            method;
            tspan=(0.0, 1.0),
            saveat=0.0,
        )

        @test_throws ArgumentError simulate(
            model,
            Up(),
            method;
            tspan=(0.0, 1.0),
            saveat=Float64[],
        )
    end
end
