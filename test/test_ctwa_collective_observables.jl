@testset "CTWA collective observables" begin

    @testset "Polarized z state" begin
        model =
            SpinModel(
                Chain(4),
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
                    cluster_size=2,
                    trajectories=8,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(100),
            )

        @test expectation(
            result,
            SpinSecondMoment(:z),
        ) == [
            1.0,
            1.0,
        ]

        @test expectation(
            result,
            SpinVariance(:z),
        ) == [
            0.0,
            0.0,
        ]
    end


    @testset "Domain wall z second moment" begin
        model =
            SpinModel(
                Chain(4),
                Field(
                    :z,
                    0.0,
                ),
            )

        result =
            simulate(
                model,
                DomainWall(),
                CTWA(
                    cluster_size=2,
                    trajectories=8,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(110),
            )

        # For |up up down down>, S_z = 0 exactly.
        @test expectation(
            result,
            SpinSecondMoment(:z),
        ) == [
            0.0,
            0.0,
        ]

        @test expectation(
            result,
            SpinVariance(:z),
        ) == [
            0.0,
            0.0,
        ]
    end


    @testset "Second moment is assembled from physical correlations" begin
        model =
            SpinModel(
                Chain(4),
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
                    cluster_size=2,
                    trajectories=6,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(120),
            )

        manual =
            fill(
                4.0,
                length(result.t),
            )

        for i in 1:3
            for j in (i + 1):4
                manual .+=
                    2.0 .*
                    expectation(
                        result,
                        Correlation(
                            :z,
                            i,
                            :z,
                            j,
                        ),
                    )
            end
        end

        manual ./=
            16.0

        @test expectation(
            result,
            SpinSecondMoment(:z),
        ) == manual
    end


    @testset "Cluster-size-one collective formula" begin
        model =
            SpinModel(
                Chain(3),
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
                    trajectories=7,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(130),
            )

        @test expectation(
            result,
            SpinSecondMoment(:z),
        ) == [
            1.0,
            1.0,
        ]

        @test expectation(
            result,
            SpinVariance(:z),
        ) == [
            0.0,
            0.0,
        ]
    end


    @testset "Float32 preservation" begin
        model =
            SpinModel(
                Chain(2),
                Field(
                    :z,
                    0.0f0,
                ),
            )

        result =
            simulate(
                model,
                Up(),
                CTWA(
                    cluster_size=2,
                    trajectories=4,
                );
                tspan=(0.0f0, 0.1f0),
                saveat=Float32[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(140),
                T=Float32,
            )

        second =
            expectation(
                result,
                SpinSecondMoment(:z),
            )

        variance =
            expectation(
                result,
                SpinVariance(:z),
            )

        @test eltype(second) === Float32
        @test eltype(variance) === Float32

        @test second == Float32[
            1.0,
            1.0,
        ]

        @test variance == Float32[
            0.0,
            0.0,
        ]
    end
end
