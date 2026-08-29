@testset "CTWA observables" begin

    @testset "Local magnetization maps physical sites to cluster generators" begin
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
                rng=Xoshiro(10),
            )

        for site in 1:4
            values =
                expectation(
                    result,
                    LocalMagnetization(
                        :z,
                        site,
                    ),
                )

            @test values == [
                1.0,
                1.0,
            ]
        end
    end


    @testset "Magnetization density" begin
        model =
            SpinModel(
                Chain(4),
                Field(
                    :z,
                    0.0,
                ),
            )

        up =
            simulate(
                model,
                Up(),
                CTWA(
                    cluster_size=2,
                    trajectories=5,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(11),
            )

        down =
            simulate(
                model,
                Down(),
                CTWA(
                    cluster_size=2,
                    trajectories=5,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(12),
            )

        @test expectation(
            up,
            Magnetization(:z),
        ) == [
            1.0,
            1.0,
        ]

        @test expectation(
            down,
            Magnetization(:z),
        ) == [
            -1.0,
            -1.0,
        ]
    end


    @testset "Same-cluster correlation uses one Pauli-string coordinate" begin
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
                    cluster_size=2,
                    trajectories=6,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(20),
            )

        values =
            expectation(
                result,
                Correlation(
                    :z,
                    1,
                    :z,
                    2,
                ),
            )

        @test values == [
            1.0,
            1.0,
        ]

        generator =
            pauli_index(
                result.basis,
                (
                    3,
                    3,
                ),
            )

        direct =
            vec(
                sum(
                    result.trajectories[
                        generator,
                        1,
                        :,
                        :,
                    ];
                    dims=2,
                ),
            ) ./
            size(
                result.trajectories,
                4,
            )

        @test values == direct
    end


    @testset "Inter-cluster correlation uses product of cluster coordinates" begin
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
                    trajectories=7,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(30),
            )

        values =
            expectation(
                result,
                Correlation(
                    :z,
                    2,
                    :z,
                    3,
                ),
            )

        @test values == [
            1.0,
            1.0,
        ]

        left =
            local_pauli_index(
                result.basis,
                result.clustering,
                2,
                3,
            )

        right =
            local_pauli_index(
                result.basis,
                result.clustering,
                3,
                3,
            )

        direct =
            zeros(
                Float64,
                length(result.t),
            )

        for trajectory in 1:size(result.trajectories, 4)
            direct .+=
                result.trajectories[
                    left,
                    1,
                    :,
                    trajectory,
                ] .*
                result.trajectories[
                    right,
                    2,
                    :,
                    trajectory,
                ]
        end

        direct ./=
            size(
                result.trajectories,
                4,
            )

        @test values == direct
    end


    @testset "Same-site Pauli identity" begin
        model =
            SpinModel(
                Chain(2),
                Field(
                    :x,
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
                rng=Xoshiro(40),
            )

        @test expectation(
            result,
            Correlation(
                :z,
                1,
                :z,
                1,
            ),
        ) == ones(3)

        @test_throws ArgumentError expectation(
            result,
            Correlation(
                :x,
                1,
                :y,
                1,
            ),
        )
    end


    @testset "Connected correlation subtracts one-point functions" begin
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
                rng=Xoshiro(50),
            )

        same_cluster =
            expectation(
                result,
                ConnectedCorrelation(
                    :z,
                    1,
                    :z,
                    2,
                ),
            )

        inter_cluster =
            expectation(
                result,
                ConnectedCorrelation(
                    :z,
                    2,
                    :z,
                    3,
                ),
            )

        @test same_cluster == zeros(2)
        @test inter_cluster == zeros(2)
    end


    @testset "Domain wall site mapping" begin
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
                    trajectories=5,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(60),
            )

        @test expectation(
            result,
            LocalMagnetization(
                :z,
                1,
            ),
        ) == [
            1.0,
            1.0,
        ]

        @test expectation(
            result,
            LocalMagnetization(
                :z,
                4,
            ),
        ) == [
            -1.0,
            -1.0,
        ]

        @test expectation(
            result,
            Magnetization(:z),
        ) == [
            0.0,
            0.0,
        ]
    end


    @testset "Invalid observable indices and axes" begin
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
                    trajectories=2,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(70),
            )

        @test_throws BoundsError expectation(
            result,
            LocalMagnetization(
                :z,
                3,
            ),
        )

        @test_throws ArgumentError Correlation(
            :z,
            0,
            :z,
            1,
        )

        @test_throws BoundsError expectation(
            result,
            Correlation(
                :z,
                3,
                :z,
                1,
            ),
        )
    end
end
