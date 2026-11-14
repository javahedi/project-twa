@testset "CTWA method and compilation" begin

    @testset "Method configuration" begin
        method =
            CTWA(
                cluster_size=3,
                trajectories=250,
            )

        @test method.cluster_size == 3
        @test method.trajectories == 250

        default =
            CTWA()

        @test default.cluster_size == 2
        @test default.trajectories == 1000

        @test_throws ArgumentError CTWA(
            cluster_size=0,
        )

        @test_throws ArgumentError CTWA(
            trajectories=0,
        )
    end


    @testset "Compiled object bundles CTWA setup" begin
        model =
            SpinModel(
                Chain(4),
                XXZ(
                    J=1.0,
                    Δ=0.5,
                ) +
                Field(
                    :z,
                    0.3,
                ),
            )

        state =
            Up()

        method =
            CTWA(
                cluster_size=2,
                trajectories=10,
            )

        compiled =
            compile(
                model,
                state,
                method,
            )

        @test compiled isa CompiledCTWA

        @test clustering(compiled).nsites == 4
        @test clustering(compiled).cluster_size == 2
        @test cluster_count(clustering(compiled)) == 2

        @test basis_size(basis(compiled)) == 15

        @test compiled.hamiltonian.basis === basis(compiled)
        @test compiled.sampling.basis.cluster_size ==
              basis(compiled).cluster_size

        @test compiled.sampling.clustering.nsites ==
              clustering(compiled).nsites
    end


    @testset "Compiled sample convenience API" begin
        model =
            SpinModel(
                Chain(4),
                Field(
                    :z,
                    0.2,
                ),
            )

        compiled =
            compile(
                model,
                Up(),
                CTWA(
                    cluster_size=2,
                    trajectories=5,
                ),
            )

        x =
            sample_ctwa(
                compiled;
                rng=Xoshiro(123),
            )

        @test size(x) ==
              (
                  basis_size(basis(compiled)),
                  cluster_count(clustering(compiled)),
              )

        y =
            similar(x)

        returned =
            sample_ctwa!(
                y,
                compiled;
                rng=Xoshiro(456),
            )

        @test returned === y
    end


    @testset "Compiled RHS delegates to cached production path" begin
        model =
            SpinModel(
                Chain(4),
                XXZ(
                    J=1.0,
                    Δ=0.5,
                ) +
                Field(
                    :z,
                    0.3,
                ),
            )

        compiled =
            compile(
                model,
                Up(),
                CTWA(
                    cluster_size=2,
                    trajectories=4,
                ),
            )

        x =
            sample_ctwa(
                compiled;
                rng=Xoshiro(7),
            )

        dx_compiled =
            similar(x)

        dx_explicit =
            similar(x)

        ctwa_rhs!(
            dx_compiled,
            x,
            compiled,
            0.0,
        )

        ctwa_rhs!(
            dx_explicit,
            x,
            compiled.hamiltonian,
            compiled.algebra,
            0.0,
        )

        @test dx_compiled == dx_explicit
    end


    @testset "Float32 compilation" begin
        model =
            SpinModel(
                Chain(2),
                XXZ(
                    J=1.0f0,
                    Δ=0.5f0,
                ),
            )

        compiled =
            compile(
                model,
                Up(),
                CTWA(
                    cluster_size=1,
                    trajectories=3,
                );
                T=Float32,
            )

        
        @test eltype(compiled.hamiltonian.local_terms) <:
              CTWALocalTerm{Float32}
        @test eltype(compiled.hamiltonian.intercluster_terms) <:
              CTWAInterClusterTerm{Float32}

        x =
            sample_ctwa(
                compiled;
                rng=Xoshiro(11),
            )

        dx =
            similar(x)

        ctwa_rhs!(
            dx,
            x,
            compiled,
            0.0f0,
        )

        @test eltype(x) == Float32
        @test eltype(dx) == Float32
    end


    @testset "Invalid cluster partition is rejected at compile time" begin
        model =
            SpinModel(
                Chain(5),
                XXZ(),
            )

        @test_throws ArgumentError compile(
            model,
            Up(),
            CTWA(
                cluster_size=2,
                trajectories=10,
            ),
        )
    end
end



@testset "Explicit CTWA clustering" begin
    lattice = SquareLattice(4, 4)

    clustering =
        Clustering(
            lattice,
            (2, 2),
        )

    method =
        CTWA(
            clustering=clustering,
            trajectories=10,
        )

    @test method.cluster_size == 4
    @test method.trajectories == 10
    @test method.clustering === clustering

    @test_throws ArgumentError CTWA(
        cluster_size=2,
        clustering=clustering,
    )

    matching =
        CTWA(
            cluster_size=4,
            clustering=clustering,
        )

    @test matching.cluster_size == 4
end


@testset "Compile with spatial square-lattice clustering" begin
    lattice =
        SquareLattice(
            4,
            4,
        )

    spatial_clustering =
        Clustering(
            lattice,
            (2, 2),
        )

    model =
        SpinModel(
            lattice,
            XXZ(
                J=1.0,
                Δ=0.5,
            ),
        )

    state = Down()

    method =
        CTWA(
            clustering=spatial_clustering,
            trajectories=2,
        )

    compiled =
        compile(
            model,
            state,
            method,
        )

    c = clustering(compiled)

    @test c === spatial_clustering
    @test c.cluster_size == 4
    @test length(c) == 4

    @test collect(cluster_sites(c, 1)) ==
          [1, 2, 5, 6]

    @test collect(cluster_sites(c, 2)) ==
          [3, 4, 7, 8]

    @test site_cluster_position(c, 5) ==
          (1, 3)

    @test site_cluster_position(c, 7) ==
          (2, 3)

    @test same_cluster(c, 1, 6)
    @test !same_cluster(c, 2, 3)
end