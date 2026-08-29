@testset "Approximation methods" begin

    @testset "DTWA construction" begin
        method = DTWA()

        @test method isa AbstractApproximation
        @test method.trajectories == 1000
        @test isconcretetype(typeof(method))
    end


    @testset "Custom trajectory count" begin
        method = DTWA(trajectories=250)

        @test method.trajectories == 250
        @test method.trajectories isa Int
    end


    @testset "Integer conversion" begin
        method = DTWA(trajectories=Int32(64))

        @test method.trajectories == 64
        @test method.trajectories isa Int
    end


    @testset "Validation" begin
        @test_throws ArgumentError DTWA(trajectories=0)
        @test_throws ArgumentError DTWA(trajectories=-1)
    end
end
