@testset "CTWA simulation" begin

    @testset "Basic simulation result" begin
        model =
            SpinModel(
                Chain(2),
                Field(
                    :z,
                    0.3,
                ),
            )

        result =
            simulate(
                model,
                Up(),
                CTWA(
                    cluster_size=1,
                    trajectories=4,
                );
                tspan=(0.0, 0.2),
                saveat=[
                    0.0,
                    0.1,
                    0.2,
                ],
                rng=Xoshiro(1234),
            )

        @test result isa CTWAResult

        @test result.t == [
            0.0,
            0.1,
            0.2,
        ]

        @test size(result.trajectories) ==
              (
                  3,
                  2,
                  3,
                  4,
              )

        @test result.basis.cluster_size == 1
        @test result.clustering.nsites == 2
        @test cluster_count(result.clustering) == 2
    end


    @testset "Clustered result layout" begin
        model =
            SpinModel(
                Chain(4),
                XXZ(
                    J=1.0,
                    Δ=0.5,
                ),
            )

        result =
            simulate(
                model,
                Up(),
                CTWA(
                    cluster_size=2,
                    trajectories=3,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(7),
            )

        @test size(result.trajectories) ==
              (
                  15,
                  2,
                  2,
                  3,
              )

        @test basis_size(result.basis) == 15
        @test result.clustering.cluster_size == 2
    end


    @testset "Initial saved state is sampled CTWA state" begin
        model =
            SpinModel(
                Chain(2),
                Field(
                    :z,
                    0.4,
                ),
            )

        result =
            simulate(
                model,
                Up(),
                CTWA(
                    cluster_size=1,
                    trajectories=5,
                );
                tspan=(0.0, 0.05),
                saveat=[
                    0.0,
                    0.05,
                ],
                rng=Xoshiro(99),
            )

        z =
            pauli_index(
                result.basis,
                (3,),
            )

        # Up() is a Z eigenstate, so the Z Weyl symbol is deterministic.
        @test all(
            result.trajectories[
                z,
                :,
                1,
                :,
            ] .== 1.0
        )
    end


    @testset "Zero Hamiltonian stays constant" begin
        model =
            SpinModel(
                Chain(2),
                Field(
                    :z,
                    0.0,
                ),
            )

        result =
            simulate(
                model,
                Up(),
                CTWA(
                    cluster_size=1,
                    trajectories=3,
                );
                tspan=(0.0, 0.2),
                saveat=[
                    0.0,
                    0.1,
                    0.2,
                ],
                rng=Xoshiro(12),
                abstol=1e-10,
                reltol=1e-10,
            )

        @test all(
            result.trajectories[
                :,
                :,
                1,
                :,
            ] .==
            result.trajectories[
                :,
                :,
                2,
                :,
            ]
        )

        @test all(
            result.trajectories[
                :,
                :,
                1,
                :,
            ] .==
            result.trajectories[
                :,
                :,
                3,
                :,
            ]
        )
    end


    @testset "Default save grid" begin
        model =
            SpinModel(
                Chain(2),
                Field(
                    :x,
                    0.2,
                ),
            )

        result =
            simulate(
                model,
                Up(),
                CTWA(
                    cluster_size=1,
                    trajectories=1,
                );
                tspan=(0.0, 1.0),
                rng=Xoshiro(5),
            )

        @test length(result.t) == 101
        @test first(result.t) == 0.0
        @test last(result.t) == 1.0
        @test size(result.trajectories, 3) == 101
    end


    @testset "Float32 simulation" begin
        model =
            SpinModel(
                Chain(2),
                Field(
                    :z,
                    0.2f0,
                ),
            )

        result =
            simulate(
                model,
                Up(),
                CTWA(
                    cluster_size=1,
                    trajectories=2,
                );
                tspan=(0.0f0, 0.1f0),
                saveat=Float32[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(18),
                T=Float32,
            )

        @test eltype(result.t) == Float32
        @test eltype(result.trajectories) == Float32
    end


    @testset "Invalid simulation arguments" begin
        model =
            SpinModel(
                Chain(2),
                Field(
                    :z,
                    0.2,
                ),
            )

        method =
            CTWA(
                cluster_size=1,
                trajectories=1,
            )

        @test_throws ArgumentError simulate(
            model,
            Up(),
            method;
            tspan=(0.0, 0.0),
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
            saveat=Float64[],
        )
    end
end
