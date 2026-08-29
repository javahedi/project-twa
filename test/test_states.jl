@testset "States" begin

    @testset "Polarized states" begin
        px = Polarized(:x)
        py = Polarized(:y; sign=-1)
        pz = Polarized(:z)

        @test px.axis === :x
        @test px.sign == 1

        @test py.axis === :y
        @test py.sign == -1

        @test pz.axis === :z
        @test pz.sign == 1
    end


    @testset "Convenience constructors" begin
        up = Up()
        down = Down()

        @test up isa Polarized
        @test up.axis === :z
        @test up.sign == 1

        @test down isa Polarized
        @test down.axis === :z
        @test down.sign == -1
    end


    @testset "Polarized directions" begin
        @test state_direction(
            Polarized(:x),
            2,
            5,
        ) == (1, 0, 0)

        @test state_direction(
            Polarized(:y; sign=-1),
            4,
            5,
        ) == (0, -1, 0)

        @test state_direction(
            Down(),
            1,
            5,
        ) == (0, 0, -1)
    end


    @testset "Domain wall" begin
        wall = DomainWall()

        @test wall.axis === :z

        @test state_direction(wall, 1, 4) == (0, 0, 1)
        @test state_direction(wall, 2, 4) == (0, 0, 1)
        @test state_direction(wall, 3, 4) == (0, 0, -1)
        @test state_direction(wall, 4, 4) == (0, 0, -1)
    end


    @testset "Odd domain wall" begin
        wall = DomainWall(axis=:x)

        @test state_direction(wall, 1, 5) == (1, 0, 0)
        @test state_direction(wall, 2, 5) == (1, 0, 0)
        @test state_direction(wall, 3, 5) == (1, 0, 0)
        @test state_direction(wall, 4, 5) == (-1, 0, 0)
        @test state_direction(wall, 5, 5) == (-1, 0, 0)
    end


    @testset "Validation" begin
        @test_throws ArgumentError Polarized(:q)
        @test_throws ArgumentError Polarized(:z; sign=0)
        @test_throws ArgumentError Polarized(:z; sign=2)
        @test_throws ArgumentError DomainWall(axis=:q)

        state = Up()

        @test_throws BoundsError state_direction(state, 0, 4)
        @test_throws BoundsError state_direction(state, 5, 4)
        @test_throws ArgumentError state_direction(state, 1, 0)
    end


    @testset "Concrete state types" begin
        @test isconcretetype(typeof(Polarized(:x)))
        @test isconcretetype(typeof(DomainWall()))
    end
end
