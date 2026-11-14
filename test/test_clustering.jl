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




@testset "Arbitrary cluster membership" begin
    # Four-by-four row-major lattice:
    #
    #  1  2   3  4
    #  5  6   7  8
    #
    #  9 10  11 12
    # 13 14  15 16
    #
    # Partition into spatial 2×2 blocks.

    sites = [
        1   3   9  11
        2   4  10  12
        5   7  13  15
        6   8  14  16
    ]

    clustering =
        Clustering(
            16,
            sites,
        )

    @test clustering.nsites == 16
    @test clustering.cluster_size == 4
    @test cluster_count(clustering) == 4

    @test collect(cluster_sites(clustering, 1)) ==
          [1, 2, 5, 6]

    @test collect(cluster_sites(clustering, 2)) ==
          [3, 4, 7, 8]

    @test site_cluster_position(clustering, 1) == (1, 1)
    @test site_cluster_position(clustering, 2) == (1, 2)
    @test site_cluster_position(clustering, 5) == (1, 3)
    @test site_cluster_position(clustering, 6) == (1, 4)

    @test site_cluster_position(clustering, 11) == (4, 1)
    @test site_cluster_position(clustering, 16) == (4, 4)

    @test same_cluster(clustering, 1, 6)
    @test !same_cluster(clustering, 1, 3)
end



@testset "Clustering partition validation" begin
    @test_throws ArgumentError Clustering(
        4,
        [
            1  3
            2  3
        ],
    )

    @test_throws ArgumentError Clustering(
        4,
        [
            1  3
            2  5
        ],
    )
end



@testset "Square-lattice block clustering" begin
    lattice =
        SquareLattice(
            4,
            4,
        )

    clustering =
        Clustering(
            lattice,
            (2, 2),
        )

    @test clustering.nsites == 16
    @test clustering.cluster_size == 4
    @test length(clustering) == 4

    @test collect(cluster_sites(clustering, 1)) ==
          [1, 2, 5, 6]

    @test collect(cluster_sites(clustering, 2)) ==
          [3, 4, 7, 8]

    @test collect(cluster_sites(clustering, 3)) ==
          [9, 10, 13, 14]

    @test collect(cluster_sites(clustering, 4)) ==
          [11, 12, 15, 16]

    @test site_cluster_position(
        clustering,
        5,
    ) == (1, 3)

    @test site_cluster_position(
        clustering,
        12,
    ) == (4, 2)

    @test same_cluster(
        clustering,
        1,
        6,
    )

    @test !same_cluster(
        clustering,
        2,
        3,
    )
end



@testset "Square-lattice block validation" begin
    lattice =
        SquareLattice(
            5,
            4,
        )

    @test_throws ArgumentError Clustering(
        lattice,
        (2, 2),
    )

    @test_throws ArgumentError Clustering(
        lattice,
        (1, 3),
    )

    @test_throws ArgumentError Clustering(
        lattice,
        (0, 2),
    )
end


@testset "Square-lattice benchmark cluster shapes" begin
    lattice = SquareLattice(8, 6)

    c2 = Clustering(lattice, (2, 1))
    @test c2.cluster_size == 2
    @test length(c2) == 24
    @test collect(cluster_sites(c2, 1)) == [1, 2]
    @test collect(cluster_sites(c2, 4)) == [7, 8]

    c4 = Clustering(lattice, (2, 2))
    @test c4.cluster_size == 4
    @test length(c4) == 12
    @test collect(cluster_sites(c4, 1)) == [1, 2, 9, 10]

    c6 = Clustering(lattice, (2, 3))
    @test c6.cluster_size == 6
    @test length(c6) == 8
    @test collect(cluster_sites(c6, 1)) ==
          [1, 2, 9, 10, 17, 18]

    c8 = Clustering(lattice, (4, 2))
    @test c8.cluster_size == 8
    @test length(c8) == 6
    @test collect(cluster_sites(c8, 1)) ==
          [1, 2, 3, 4, 9, 10, 11, 12]
end