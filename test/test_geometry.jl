@testset "Geometry" begin
    chain = Chain(50)

    @test chain.nsites == 50
    @test chain.boundary isa OpenBoundary

    periodic_chain = Chain(20, PeriodicBoundary())

    @test periodic_chain.nsites == 20
    @test periodic_chain.boundary isa PeriodicBoundary

    @test_throws ArgumentError Chain(0)
    @test_throws ArgumentError Chain(-1)
end

@testset "Open-chain bonds" begin
    chain = Chain(4)

    @test collect(bonds(chain)) == [
        (1, 2),
        (2, 3),
        (3, 4),
    ]

    @test length(bonds(chain)) == 3
    @test eltype(bonds(chain)) == Tuple{Int,Int}
end


@testset "Periodic-chain bonds" begin
    chain = Chain(4, PeriodicBoundary())

    @test collect(bonds(chain)) == [
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 1),
    ]

    @test length(bonds(chain)) == 4
end


@testset "Small-chain bonds" begin
    @test isempty(bonds(Chain(1)))
    @test length(bonds(Chain(1))) == 0

    @test collect(bonds(Chain(2))) == [
        (1, 2),
    ]

    @test collect(
        bonds(Chain(2, PeriodicBoundary())),
    ) == [
        (1, 2),
    ]
end


@testset "Chain distance" begin

    @testset "Open boundary" begin
        chain = Chain(6)

        @test distance(chain, 1, 1) == 0
        @test distance(chain, 1, 2) == 1
        @test distance(chain, 1, 5) == 4
        @test distance(chain, 1, 6) == 5

        # Distance should be symmetric.
        @test distance(chain, 5, 1) == 4
    end


    @testset "Periodic boundary" begin
        chain = Chain(6, PeriodicBoundary())

        @test distance(chain, 1, 1) == 0
        @test distance(chain, 1, 2) == 1

        # Going through the periodic boundary is shorter.
        @test distance(chain, 1, 6) == 1
        @test distance(chain, 1, 5) == 2

        # For an even chain, opposite sites are N/2 apart.
        @test distance(chain, 1, 4) == 3

        @test distance(chain, 5, 1) == 2
    end


    @testset "Invalid sites" begin
        chain = Chain(6)

        @test_throws BoundsError distance(chain, 0, 1)
        @test_throws BoundsError distance(chain, 1, 0)

        @test_throws BoundsError distance(chain, 7, 1)
        @test_throws BoundsError distance(chain, 1, 7)
    end
end



@testset "Geometry site-count interface" begin
    @test nsites(Chain(1)) == 1
    @test nsites(Chain(10)) == 10
    @test nsites(Chain(10, PeriodicBoundary())) == 10
end


@testset "Square lattice geometry" begin
    lattice = SquareLattice(4, 3)

    @test nsites(lattice) == 12
    @test lattice.Lx == 4
    @test lattice.Ly == 3

    @test coordinates(lattice, 1) == (1, 1)
    @test coordinates(lattice, 4) == (4, 1)
    @test coordinates(lattice, 5) == (1, 2)
    @test coordinates(lattice, 8) == (4, 2)
    @test coordinates(lattice, 12) == (4, 3)

    periodic = SquareLattice(
        4,
        3,
        PeriodicBoundary(),
    )

    @test nsites(periodic) == 12
    @test coordinates(periodic, 12) == (4, 3)

    @test_throws BoundsError coordinates(lattice, 0)
    @test_throws BoundsError coordinates(lattice, 13)
end



@testset "Square lattice bonds" begin
    open = SquareLattice(3, 2)

    @test Set(bonds(open)) == Set([
        (1, 2),
        (2, 3),
        (4, 5),
        (5, 6),
        (1, 4),
        (2, 5),
        (3, 6),
    ])

    @test length(bonds(open)) == 7

    periodic = SquareLattice(
        3,
        2,
        PeriodicBoundary(),
    )

    @test length(bonds(periodic)) == 9

    @test Set(bonds(periodic)) == Set([
        (1, 2),
        (2, 3),
        (3, 1),
        (4, 5),
        (5, 6),
        (6, 4),
        (1, 4),
        (2, 5),
        (3, 6),

    ])
end


@testset "Small periodic square lattices" begin
    lattice = SquareLattice(
        2,
        2,
        PeriodicBoundary(),
    )

    @test length(bonds(lattice)) == 4

    @test Set(bonds(lattice)) == Set([
        (1, 2),
        (3, 4),
        (1, 3),
        (2, 4),
    ])
end


@testset "Square lattice distances" begin
    lattice = SquareLattice(4, 3)

    @test distance(lattice, 1, 1) == 0
    @test distance(lattice, 1, 2) == 1
    @test distance(lattice, 1, 5) == 1
    @test distance(lattice, 1, 6) == sqrt(2)
    @test distance(lattice, 1, 12) == sqrt(13)

    periodic = SquareLattice(
        4,
        3,
        PeriodicBoundary(),
    )

    @test distance(periodic, 1, 4) == 1
    @test distance(periodic, 1, 9) == 1
    @test distance(periodic, 1, 12) == sqrt(2)

    @test_throws BoundsError distance(lattice, 0, 1)
    @test_throws BoundsError distance(lattice, 1, 13)
end