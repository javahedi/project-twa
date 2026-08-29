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