@testset "CTWA correlation profiles and structure factors" begin

    @testset "Open-chain z correlation profile for polarized state" begin
        geometry =
            Chain(4)

        model =
            SpinModel(
                geometry,
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
                rng=Xoshiro(200),
            )

        profile =
            correlation_profile(
                result,
                geometry,
                :z,
            )

        @test profile.distances ==
              collect(0:3)

        @test profile.values ==
              ones(
                  2,
                  4,
              )
    end


    @testset "Periodic profile distance range" begin
        geometry =
            Chain(
                6,
                PeriodicBoundary(),
            )

        model =
            SpinModel(
                geometry,
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
                    trajectories=4,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(210),
            )

        profile =
            correlation_profile(
                result,
                geometry,
                :z,
            )

        @test profile.distances ==
              collect(0:3)

        @test profile.values ==
              ones(
                  2,
                  4,
              )
    end


    @testset "Connected profile vanishes for product z eigenstate" begin
        geometry =
            Chain(4)

        model =
            SpinModel(
                geometry,
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
                rng=Xoshiro(220),
            )

        profile =
            correlation_profile(
                result,
                geometry,
                :z;
                connected=true,
            )

        @test profile.values ==
              zeros(
                  2,
                  4,
              )
    end


    @testset "q=0 structure factor agrees with collective second moment" begin
        geometry =
            Chain(4)

        model =
            SpinModel(
                geometry,
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
                    trajectories=7,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(230),
            )

        sf =
            structure_factor(
                result,
                geometry,
                :z,
            )

        second =
            expectation(
                result,
                SpinSecondMoment(:z),
            )

        @test sf.values[:, 1] ≈
              geometry.nsites .*
              second
    end


    @testset "Connected q=0 agrees with collective variance" begin
        geometry =
            Chain(4)

        model =
            SpinModel(
                geometry,
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
                    cluster_size=2,
                    trajectories=16,
                );
                tspan=(0.0, 0.2),
                saveat=[
                    0.0,
                    0.1,
                    0.2,
                ],
                rng=Xoshiro(240),
            )

        sf =
            structure_factor(
                result,
                geometry,
                :z;
                connected=true,
            )

        variance =
            expectation(
                result,
                SpinVariance(:z),
            )

        @test sf.values[:, 1] ≈
              geometry.nsites .*
              variance
    end


    @testset "Polarized state structure factor" begin
        geometry =
            Chain(4)

        model =
            SpinModel(
                geometry,
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
                    trajectories=5,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(250),
            )

        sf =
            structure_factor(
                result,
                geometry,
                :z,
            )

        @test sf.values[:, 1] ==
              fill(
                  4.0,
                  2,
              )

       @test isapprox(
                sf.values[:, 2:end],
                zeros(2, 3);
                atol=1e-14,
            )
    end


    @testset "Custom momenta" begin
        geometry =
            Chain(4)

        model =
            SpinModel(
                geometry,
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
                    trajectories=4,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(260),
            )

        qs =
            [
                0.0,
                pi,
            ]

        sf =
            structure_factor(
                result,
                geometry,
                :z;
                momenta=qs,
            )

        @test sf.momenta ==
              qs

        @test size(
            sf.values,
        ) == (
            2,
            2,
        )

        @test sf.values[:, 1] ==
              fill(
                  4.0,
                  2,
              )

        # @test sf.values[:, 2] ≈
        #       zeros(
        #           2,
        #       )

        @test isapprox(
                sf.values[:, 2],
                zeros(2);
                atol=1e-14,
            )
    end


    @testset "Geometry/result mismatch" begin
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
                    trajectories=2,
                );
                tspan=(0.0, 0.1),
                saveat=[
                    0.0,
                    0.1,
                ],
                rng=Xoshiro(270),
            )

        @test_throws ArgumentError correlation_profile(
            result,
            Chain(6),
            :z,
        )

        @test_throws ArgumentError structure_factor(
            result,
            Chain(6),
            :z,
        )
    end
end
