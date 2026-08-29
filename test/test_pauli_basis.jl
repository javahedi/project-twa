@testset "Pauli-string basis" begin

    @testset "Basis dimensions" begin
        basis1 = PauliStringBasis(1)
        basis2 = PauliStringBasis(2)
        basis3 = PauliStringBasis(3)

        @test basis1.cluster_size == 1
        @test basis_size(basis1) == 3
        @test basis_size(basis2) == 15
        @test basis_size(basis3) == 63

        @test length(basis1) == 3
        @test length(basis2) == 15
        @test length(basis3) == 63

        @test_throws ArgumentError PauliStringBasis(0)
        @test_throws ArgumentError PauliStringBasis(-1)
    end


    @testset "One-site ordering" begin
        basis = PauliStringBasis(1)

        @test pauli_digits(basis, 1) == (1,)
        @test pauli_digits(basis, 2) == (2,)
        @test pauli_digits(basis, 3) == (3,)

        @test pauli_symbols(basis, 1) == (:X,)
        @test pauli_symbols(basis, 2) == (:Y,)
        @test pauli_symbols(basis, 3) == (:Z,)

        @test pauli_index(basis, (1,)) == 1
        @test pauli_index(basis, (2,)) == 2
        @test pauli_index(basis, (3,)) == 3
    end


    @testset "Two-site canonical ordering" begin
        basis = PauliStringBasis(2)

        # Base-4 lexicographic ordering:
        #
        # 1 = 01 = IX
        # 2 = 02 = IY
        # 3 = 03 = IZ
        # 4 = 10 = XI
        # 5 = 11 = XX
        @test pauli_digits(basis, 1) == (0, 1)
        @test pauli_digits(basis, 2) == (0, 2)
        @test pauli_digits(basis, 3) == (0, 3)
        @test pauli_digits(basis, 4) == (1, 0)
        @test pauli_digits(basis, 5) == (1, 1)

        @test pauli_symbols(basis, 1) == (:I, :X)
        @test pauli_symbols(basis, 4) == (:X, :I)
        @test pauli_symbols(basis, 5) == (:X, :X)

        @test pauli_index(basis, (0, 1)) == 1
        @test pauli_index(basis, (1, 0)) == 4
        @test pauli_index(basis, (3, 3)) == 15
    end


    @testset "Code and index distinction" begin
        basis = PauliStringBasis(2)

        # The identity has a valid raw Pauli code but is not a traceless
        # basis element.
        @test pauli_code(basis, (0, 0)) == 0

        @test_throws ArgumentError pauli_index(
            basis,
            (0, 0),
        )

        @test pauli_code(basis, (2, 3)) == 11
        @test pauli_index(basis, (2, 3)) == 11
    end


    @testset "Round-trip mapping" begin
        for n in 1:4
            basis = PauliStringBasis(n)

            @test all(
                pauli_index(
                    basis,
                    pauli_digits(basis, index),
                ) == index
                for index in 1:basis_size(basis)
            )
        end
    end


    @testset "Three-site examples" begin
        basis = PauliStringBasis(3)

        # IIX = 0*16 + 0*4 + 1
        @test pauli_index(
            basis,
            (0, 0, 1),
        ) == 1

        # XII = 1*16
        @test pauli_index(
            basis,
            (1, 0, 0),
        ) == 16

        # XYZ = 1*16 + 2*4 + 3
        @test pauli_index(
            basis,
            (1, 2, 3),
        ) == 27

        @test pauli_digits(
            basis,
            27,
        ) == (1, 2, 3)
    end


    @testset "Input validation" begin
        basis = PauliStringBasis(2)

        @test_throws DimensionMismatch pauli_code(
            basis,
            (1,),
        )

        @test_throws DimensionMismatch pauli_code(
            basis,
            (1, 2, 3),
        )

        @test_throws ArgumentError pauli_code(
            basis,
            (1, 4),
        )

        @test_throws ArgumentError pauli_code(
            basis,
            (1, -1),
        )

        @test_throws ArgumentError pauli_code(
            basis,
            (1, 2.0),
        )

        @test_throws BoundsError pauli_digits(
            basis,
            0,
        )

        @test_throws BoundsError pauli_digits(
            basis,
            16,
        )
    end


    @testset "Pauli symbols" begin
        @test pauli_symbol(0) == :I
        @test pauli_symbol(1) == :X
        @test pauli_symbol(2) == :Y
        @test pauli_symbol(3) == :Z

        @test_throws ArgumentError pauli_symbol(4)
        @test_throws ArgumentError pauli_symbol(-1)
    end
end
