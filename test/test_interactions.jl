@testset "Interactions" begin

    @testset "XXZ" begin
        interaction = XXZ(J=1.0, Δ=0.5)

        @test interaction.J == 1.0
        @test interaction.Δ == 0.5

        isotropic = XXZ(J=2.0)

        @test isotropic.J == 2.0
        @test isotropic.Δ == 1.0
    end


    @testset "Nearest neighbor" begin
        interaction = XXZ(J=1.0, Δ=0.5)
        coupling = NearestNeighbor(interaction)

        @test coupling.interaction === interaction
    end


    @testset "Power law" begin
        interaction = XXZ(J=1.0, Δ=0.5)

        coupling = PowerLaw(
            interaction;
            α=3.0,
        )

        @test coupling.interaction === interaction
        @test coupling.α == 3.0

        @test_throws ArgumentError PowerLaw(
            interaction;
            α=-1.0,
        )
    end


    @testset "Numeric promotion" begin
        interaction = XXZ(
            J=1,
            Δ=0.5,
        )

        @test interaction.J isa Float64
        @test interaction.Δ isa Float64

        @test interaction.J == 1.0
        @test interaction.Δ == 0.5
    end
end
