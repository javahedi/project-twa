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
