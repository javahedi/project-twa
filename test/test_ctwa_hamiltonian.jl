@testset "CTWA Hamiltonian compilation" begin

    function local_term_dict(compiled)
        Dict(
            (
                term.cluster,
                term.generator,
            ) => term.coefficient
            for term in local_terms(compiled)
        )
    end

    function inter_term_dict(compiled)
        Dict(
            (
                term.left_cluster,
                term.left_generator,
                term.right_cluster,
                term.right_generator,
            ) => term.coefficient
            for term in intercluster_terms(compiled)
        )
    end


    @testset "Nearest-neighbor XXZ with two-site clusters" begin
        model = SpinModel(
            Chain(4),
            XXZ(
                J=2.0,
                Δ=0.5,
            ),
        )

        clustering = Clustering(
            4,
            2,
        )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )

        basis = compiled.basis

        xx = pauli_index(
            basis,
            (1, 1),
        )

        yy = pauli_index(
            basis,
            (2, 2),
        )

        zz = pauli_index(
            basis,
            (3, 3),
        )

        ix = pauli_index(
            basis,
            (0, 1),
        )

        iy = pauli_index(
            basis,
            (0, 2),
        )

        iz = pauli_index(
            basis,
            (0, 3),
        )

        xi = pauli_index(
            basis,
            (1, 0),
        )

        yi = pauli_index(
            basis,
            (2, 0),
        )

        zi = pauli_index(
            basis,
            (3, 0),
        )

        local_ = local_term_dict(compiled)

        inter = inter_term_dict(
            compiled,
        )

        # Bonds (1,2) and (3,4) are intra-cluster. Therefore each XXZ
        # interaction becomes a linear term in the corresponding cluster.
        @test length(local_) == 6

        @test local_[(1, xx)] == 2.0
        @test local_[(1, yy)] == 2.0
        @test local_[(1, zz)] == 1.0

        @test local_[(2, xx)] == 2.0
        @test local_[(2, yy)] == 2.0
        @test local_[(2, zz)] == 1.0

        # Bond (2,3) crosses the cluster boundary:
        #
        # site 2 -> IX in cluster 1
        # site 3 -> XI in cluster 2
        @test length(inter) == 3

        @test inter[
            (
                1,
                ix,
                2,
                xi,
            )
        ] == 2.0

        @test inter[
            (
                1,
                iy,
                2,
                yi,
            )
        ] == 2.0

        @test inter[
            (
                1,
                iz,
                2,
                zi,
            )
        ] == 1.0
    end


    @testset "Uniform field becomes local one-body generators" begin
        model = SpinModel(
            Chain(4),
            Field(
                :z,
                0.3,
            ),
        )

        clustering = Clustering(
            4,
            2,
        )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )

        basis = compiled.basis

        zi = pauli_index(
            basis,
            (3, 0),
        )

        iz = pauli_index(
            basis,
            (0, 3),
        )

        local_ = local_term_dict(
            compiled,
        )

        @test isempty(
            intercluster_terms(compiled),
        )

        @test length(local_) == 4

        @test local_[(1, zi)] == 0.3
        @test local_[(1, iz)] == 0.3
        @test local_[(2, zi)] == 0.3
        @test local_[(2, iz)] == 0.3
    end


    @testset "Interaction plus field" begin
        model = SpinModel(
            Chain(4),
            XXZ(
                J=1.0,
                Δ=2.0,
            ) +
            Field(
                :x,
                0.25,
            ),
        )

        clustering = Clustering(
            4,
            2,
        )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )

        basis = compiled.basis
        local_ = local_term_dict(
            compiled,
        )

        xi = pauli_index(
            basis,
            (1, 0),
        )

        ix = pauli_index(
            basis,
            (0, 1),
        )

        xx = pauli_index(
            basis,
            (1, 1),
        )

        zz = pauli_index(
            basis,
            (3, 3),
        )

        @test local_[(1, xi)] == 0.25
        @test local_[(1, ix)] == 0.25
        @test local_[(2, xi)] == 0.25
        @test local_[(2, ix)] == 0.25

        @test local_[(1, xx)] == 1.0
        @test local_[(2, xx)] == 1.0

        @test local_[(1, zz)] == 2.0
        @test local_[(2, zz)] == 2.0
    end


    @testset "Single-site clusters recover pairwise structure" begin
        model = SpinModel(
            Chain(3),
            XXZ(
                J=1.5,
                Δ=0.4,
            ),
        )

        clustering = Clustering(
            3,
            1,
        )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )

        # With cluster size one, no physical pair can be intra-cluster.
        @test isempty(
            local_terms(compiled),
        )

        # Open three-site chain has two bonds, with three XXZ components
        # on each bond.
        @test length(
            intercluster_terms(compiled),
        ) == 6

        basis = compiled.basis

        x = pauli_index(
            basis,
            (1,),
        )

        y = pauli_index(
            basis,
            (2,),
        )

        z = pauli_index(
            basis,
            (3,),
        )

        inter = inter_term_dict(
            compiled,
        )

        for pair in (
            (1, 2),
            (2, 3),
        )
            left, right = pair

            @test inter[
                (
                    left,
                    x,
                    right,
                    x,
                )
            ] == 1.5

            @test inter[
                (
                    left,
                    y,
                    right,
                    y,
                )
            ] == 1.5

            @test inter[
                (
                    left,
                    z,
                    right,
                    z,
                )
            ] ≈ 0.6
        end
    end


    @testset "Periodic boundary crossing" begin
        model = SpinModel(
            Chain(
                4,
                PeriodicBoundary(),
            ),
            XXZ(
                J=1.0,
                Δ=1.0,
            ),
        )

        clustering = Clustering(
            4,
            2,
        )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )

        basis = compiled.basis
        inter = inter_term_dict(
            compiled,
        )

        # Internal nearest-neighbor bonds:
        # (1,2), (3,4)
        #
        # Cross-cluster bonds:
        # (2,3), (1,4)
        #
        # Therefore there are two cluster-boundary bonds, each carrying
        # X, Y, Z components.
        @test length(inter) == 6

        xi = pauli_index(
            basis,
            (1, 0),
        )

        ix = pauli_index(
            basis,
            (0, 1),
        )

        # Periodic bond (1,4):
        #
        # site 1 -> XI in cluster 1
        # site 4 -> IX in cluster 2
        @test inter[
            (
                1,
                xi,
                2,
                ix,
            )
        ] == 1.0
    end


    @testset "Power-law interaction creates intra and inter terms" begin
        model = SpinModel(
            Chain(4),
            PowerLaw(
                XXZ(
                    J=1.0,
                    Δ=1.0,
                );
                α=2.0,
            ),
        )

        clustering = Clustering(
            4,
            2,
        )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )

        basis = compiled.basis
        local_ = local_term_dict(
            compiled,
        )
        inter = inter_term_dict(
            compiled,
        )

        xx = pauli_index(
            basis,
            (1, 1),
        )

        # Sites 1 and 2 are distance one and belong to cluster 1.
        @test local_[(1, xx)] == 1.0

        xi = pauli_index(
            basis,
            (1, 0),
        )

        # Sites 1 and 3 are distance two:
        #
        # site 1 -> XI in cluster 1
        # site 3 -> XI in cluster 2
        @test inter[
            (
                1,
                xi,
                2,
                xi,
            )
        ] == 1 / 4
    end


    @testset "Coefficient precision" begin
        model = SpinModel(
            Chain(4),
            XXZ(
                J=Float32(1),
                Δ=Float32(0.5),
            ),
        )

        clustering = Clustering(
            4,
            2,
        )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering;
                T=Float32,
            )

        @test eltype(
            compiled.local_terms,
        ) == CTWALocalTerm{Float32}

        @test eltype(
            compiled.intercluster_terms,
        ) == CTWAInterClusterTerm{Float32}

        @test all(
            term.coefficient isa Float32
            for term in compiled.local_terms
        )

        @test all(
            term.coefficient isa Float32
            for term in compiled.intercluster_terms
        )
    end


    @testset "Model/clustering size mismatch" begin
        model = SpinModel(
            Chain(4),
            XXZ(),
        )

        clustering = Clustering(
            6,
            2,
        )

        @test_throws DimensionMismatch compile_ctwa_hamiltonian(
            model,
            clustering,
        )
    end
end
