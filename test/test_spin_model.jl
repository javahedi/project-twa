@testset "SpinModel" begin

    @testset "Bare interaction becomes nearest-neighbor" begin
        chain = Chain(50)
        xxz = XXZ(J=1.0, Δ=0.5)

        model = SpinModel(chain, xxz)

        @test geometry(model) === chain
        @test hamiltonian(model) === model.hamiltonian
        @test length(terms(model)) == 1

        pairterm = terms(model)[1]

        @test pairterm isa NearestNeighbor
        @test pairterm.interaction === xxz
        @test pairterm.interaction.J == 1.0
        @test pairterm.interaction.Δ == 0.5
    end


    @testset "Explicit nearest-neighbor profile" begin
        chain = Chain(20)
        profile = NearestNeighbor(
            XXZ(J=2.0, Δ=0.25),
        )

        model = SpinModel(chain, profile)

        @test geometry(model) === chain
        @test length(terms(model)) == 1
        @test terms(model)[1] === profile
    end


    @testset "Power-law model" begin
        chain = Chain(30, PeriodicBoundary())

        profile = PowerLaw(
            XXZ(J=1.5, Δ=0.75);
            α=3.0,
        )

        model = SpinModel(chain, profile)

        @test geometry(model) === chain
        @test length(terms(model)) == 1
        @test terms(model)[1] === profile
        @test terms(model)[1].α == 3.0
        @test terms(model)[1].interaction.J == 1.5
        @test terms(model)[1].interaction.Δ == 0.75
    end


    @testset "Interaction plus field" begin
        chain = Chain(12)
        xxz = XXZ(J=1.0, Δ=0.5)
        field = Field(:z, 0.3)

        model = SpinModel(
            chain,
            xxz + field,
        )

        model_terms = terms(model)

        @test length(model_terms) == 2
        @test model_terms[1] isa NearestNeighbor
        @test model_terms[1].interaction === xxz
        @test model_terms[2] === field
    end


    @testset "Power law plus field" begin
        chain = Chain(12, PeriodicBoundary())

        profile = PowerLaw(
            XXZ(J=1.0, Δ=0.5);
            α=2.5,
        )

        field = Field(:x, 0.2)

        model = SpinModel(
            chain,
            profile + field,
        )

        model_terms = terms(model)

        @test length(model_terms) == 2
        @test model_terms[1] === profile
        @test model_terms[2] === field
    end


    @testset "Multiple fields" begin
        xxz = XXZ()
        fx = Field(:x, 0.1)
        fz = Field(:z, 0.2)

        model = SpinModel(
            Chain(8),
            xxz + fx + fz,
        )

        model_terms = terms(model)

        @test length(model_terms) == 3
        @test model_terms[1] isa NearestNeighbor
        @test model_terms[2] === fx
        @test model_terms[3] === fz
    end


    @testset "Concrete storage" begin
        model = SpinModel(
            Chain(10),
            XXZ(J=1.0, Δ=0.5) + Field(:z, 0.2),
        )

        @test isconcretetype(typeof(model))
        @test model.hamiltonian.terms isa Tuple
    end
end
