@testset "CTWA clustering" begin

    @testset "Construction" begin
        clustering = Clustering(
            6,
            2,
        )

        @test clustering.nsites == 6
        @test clustering.cluster_size == 2
        @test clustering.nclusters == 3

        @test cluster_count(
            clustering,
        ) == 3

        @test length(
            clustering,
        ) == 3

        @test Clustering(
            Chain(6),
            2,
        ).nclusters == 3

        @test_throws ArgumentError Clustering(
            0,
            1,
        )

        @test_throws ArgumentError Clustering(
            4,
            0,
        )

        @test_throws ArgumentError Clustering(
            4,
            5,
        )

        @test_throws ArgumentError Clustering(
            5,
            2,
        )
    end


    @testset "Cluster site ranges" begin
        clustering = Clustering(
            6,
            2,
        )

        @test cluster_sites(
            clustering,
            1,
        ) == 1:2

        @test cluster_sites(
            clustering,
            2,
        ) == 3:4

        @test cluster_sites(
            clustering,
            3,
        ) == 5:6

        @test_throws BoundsError cluster_sites(
            clustering,
            0,
        )

        @test_throws BoundsError cluster_sites(
            clustering,
            4,
        )
    end


    @testset "Physical site mapping" begin
        clustering = Clustering(
            6,
            2,
        )

        expected_clusters = (
            1,
            1,
            2,
            2,
            3,
            3,
        )

        expected_positions = (
            1,
            2,
            1,
            2,
            1,
            2,
        )

        @test all(
            site_cluster(
                clustering,
                site,
            ) == expected_clusters[site]
            for site in 1:6
        )

        @test all(
            site_position(
                clustering,
                site,
            ) == expected_positions[site]
            for site in 1:6
        )

        @test all(
            site_cluster_position(
                clustering,
                site,
            ) == (
                expected_clusters[site],
                expected_positions[site],
            )
            for site in 1:6
        )

        @test_throws BoundsError site_cluster(
            clustering,
            0,
        )

        @test_throws BoundsError site_position(
            clustering,
            7,
        )
    end


    @testset "Same-cluster checks" begin
        clustering = Clustering(
            8,
            2,
        )

        @test same_cluster(
            clustering,
            1,
            2,
        )

        @test same_cluster(
            clustering,
            3,
            4,
        )

        @test !same_cluster(
            clustering,
            2,
            3,
        )

        @test !same_cluster(
            clustering,
            1,
            8,
        )
    end


    @testset "Single-site clusters" begin
        clustering = Clustering(
            4,
            1,
        )

        @test clustering.nclusters == 4

        @test all(
            cluster_sites(
                clustering,
                cluster,
            ) == cluster:cluster
            for cluster in 1:4
        )

        @test all(
            site_cluster(
                clustering,
                site,
            ) == site
            for site in 1:4
        )

        @test all(
            site_position(
                clustering,
                site,
            ) == 1
            for site in 1:4
        )
    end


    @testset "One-site physical operator mapping" begin
        clustering = Clustering(
            6,
            2,
        )

        basis = PauliStringBasis(2)

        # First site in a cluster:
        #
        # X -> XI
        @test local_pauli_digits(
            clustering,
            1,
            1,
        ) == (1, 0)

        # Second site in a cluster:
        #
        # X -> IX
        @test local_pauli_digits(
            clustering,
            2,
            1,
        ) == (0, 1)

        # The same local pattern repeats in later clusters.
        @test local_pauli_digits(
            clustering,
            3,
            2,
        ) == (2, 0)

        @test local_pauli_digits(
            clustering,
            4,
            3,
        ) == (0, 3)

        @test local_pauli_index(
            basis,
            clustering,
            1,
            1,
        ) == pauli_index(
            basis,
            (1, 0),
        )

        @test local_pauli_index(
            basis,
            clustering,
            2,
            1,
        ) == pauli_index(
            basis,
            (0, 1),
        )

        @test local_pauli_index(
            basis,
            clustering,
            4,
            3,
        ) == pauli_index(
            basis,
            (0, 3),
        )
    end


    @testset "Three-site cluster mapping" begin
        clustering = Clustering(
            6,
            3,
        )

        basis = PauliStringBasis(3)

        @test local_pauli_digits(
            clustering,
            1,
            1,
        ) == (1, 0, 0)

        @test local_pauli_digits(
            clustering,
            2,
            2,
        ) == (0, 2, 0)

        @test local_pauli_digits(
            clustering,
            3,
            3,
        ) == (0, 0, 3)

        @test local_pauli_index(
            basis,
            clustering,
            3,
            3,
        ) == pauli_index(
            basis,
            (0, 0, 3),
        )
    end


    @testset "Basis/clustering compatibility" begin
        clustering = Clustering(
            6,
            2,
        )

        wrong_basis = PauliStringBasis(3)

        @test_throws DimensionMismatch local_pauli_index(
            wrong_basis,
            clustering,
            1,
            1,
        )
    end


    @testset "Operator validation" begin
        clustering = Clustering(
            4,
            2,
        )

        @test_throws ArgumentError local_pauli_digits(
            clustering,
            1,
            0,
        )

        @test_throws ArgumentError local_pauli_digits(
            clustering,
            1,
            4,
        )

        @test_throws BoundsError local_pauli_digits(
            clustering,
            5,
            1,
        )
    end
end
