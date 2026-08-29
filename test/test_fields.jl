@testset "Fields" begin

    @testset "Construction" begin
        fx = Field(:x, 1.0)
        fy = Field(:y, 2.0)
        fz = Field(:z, 0.3)

        @test fx.axis === :x
        @test fy.axis === :y
        @test fz.axis === :z

        @test fx.strength == 1.0
        @test fy.strength == 2.0
        @test fz.strength == 0.3
    end


    @testset "Cartesian components" begin
        @test field_components(Field(:x, 2.0)) == (
            2.0,
            0.0,
            0.0,
        )

        @test field_components(Field(:y, 2.0)) == (
            0.0,
            2.0,
            0.0,
        )

        @test field_components(Field(:z, 2.0)) == (
            0.0,
            0.0,
            2.0,
        )
    end


    @testset "Invalid axis" begin
        @test_throws ArgumentError Field(:a, 1.0)
        @test_throws ArgumentError Field(:xy, 1.0)
    end


    @testset "Numeric type preservation" begin
        field = Field(:z, Float32(0.5))

        @test field.strength isa Float32

        components = field_components(field)

        @test components isa NTuple{3,Float32}
        @test components == (
            0.0f0,
            0.0f0,
            0.5f0,
        )
    end


    @testset "Integer field" begin
        field = Field(:x, 2)

        @test field.strength isa Int
        @test field_components(field) == (2, 0, 0)
    end
end
