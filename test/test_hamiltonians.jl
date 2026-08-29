@testset "Hamiltonian composition" begin

    @testset "Interaction plus field" begin
        xxz = XXZ(J=1.0, Δ=0.5)
        field = Field(:z, 0.3)

        h = xxz + field

        @test h isa Hamiltonian
        @test length(h) == 2
        @test h.terms[1] === xxz
        @test h.terms[2] === field
    end


    @testset "Coupling profile plus field" begin
        profile = PowerLaw(
            XXZ(J=1.0, Δ=0.5);
            α=3.0,
        )
        field = Field(:z, 0.3)

        h = profile + field

        @test h isa Hamiltonian
        @test length(h) == 2
        @test h.terms[1] === profile
        @test h.terms[2] === field
    end


    @testset "Multiple terms" begin
        xxz = XXZ(J=1.0, Δ=0.5)
        fx = Field(:x, 0.2)
        fz = Field(:z, 0.3)

        h = xxz + fx + fz

        @test length(h) == 3
        @test h.terms == (xxz, fx, fz)
    end


    @testset "Hamiltonian plus Hamiltonian" begin
        h1 = XXZ() + Field(:x, 0.2)
        h2 = Field(:y, 0.4) + Field(:z, 0.3)

        h = h1 + h2

        @test length(h) == 4
        @test h.terms == (
            h1.terms[1],
            h1.terms[2],
            h2.terms[1],
            h2.terms[2],
        )
    end


    @testset "Iteration" begin
        xxz = XXZ()
        fx = Field(:x, 0.1)
        fz = Field(:z, 0.2)

        h = xxz + fx + fz

        @test collect(h) == [xxz, fx, fz]
    end


    @testset "Tuple storage" begin
        h = XXZ() + Field(:z, 0.3)

        @test h.terms isa Tuple
        @test isconcretetype(typeof(h))
    end
end
