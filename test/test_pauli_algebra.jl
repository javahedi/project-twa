@testset "Pauli-string algebra" begin

    @testset "Single-site multiplication" begin
        @test pauli_product_digit(0, 0) == (1, 0)
        @test pauli_product_digit(0, 1) == (1, 1)
        @test pauli_product_digit(2, 0) == (1, 2)

        @test pauli_product_digit(1, 1) == (1, 0)
        @test pauli_product_digit(2, 2) == (1, 0)
        @test pauli_product_digit(3, 3) == (1, 0)

        @test pauli_product_digit(1, 2) == (im, 3)
        @test pauli_product_digit(2, 3) == (im, 1)
        @test pauli_product_digit(3, 1) == (im, 2)

        @test pauli_product_digit(2, 1) == (-im, 3)
        @test pauli_product_digit(3, 2) == (-im, 1)
        @test pauli_product_digit(1, 3) == (-im, 2)

        @test_throws ArgumentError pauli_product_digit(4, 1)
        @test_throws ArgumentError pauli_product_digit(1, -1)
    end


    @testset "One-site string products" begin
        basis = PauliStringBasis(1)

        x = pauli_index(basis, (1,))
        y = pauli_index(basis, (2,))
        z = pauli_index(basis, (3,))

        @test pauli_product(
            basis,
            x,
            y,
        ) == (im, z)

        @test pauli_product(
            basis,
            y,
            x,
        ) == (-im, z)

        @test pauli_product(
            basis,
            x,
            x,
        ) == (1 + 0im, 0)
    end


    @testset "Two-site products" begin
        basis = PauliStringBasis(2)

        xi = pauli_index(
            basis,
            (1, 0),
        )

        yz = pauli_index(
            basis,
            (2, 3),
        )

        zz = pauli_code(
            basis,
            (3, 3),
        )

        phase, code = pauli_product(
            basis,
            xi,
            yz,
        )

        @test phase == im
        @test code == zz

        ix = pauli_index(
            basis,
            (0, 1),
        )

        @test pauli_product(
            basis,
            xi,
            ix,
        ) == (
            1 + 0im,
            pauli_code(basis, (1, 1)),
        )

        @test pauli_product(
            basis,
            xi,
            xi,
        ) == (
            1 + 0im,
            0,
        )
    end


    @testset "Commutators" begin
        basis1 = PauliStringBasis(1)

        x = pauli_index(basis1, (1,))
        y = pauli_index(basis1, (2,))
        z = pauli_index(basis1, (3,))

        @test pauli_commutator(
            basis1,
            x,
            y,
        ) == (2im, z)

        @test pauli_commutator(
            basis1,
            y,
            x,
        ) == (-2im, z)

        @test pauli_commutator(
            basis1,
            x,
            x,
        ) == (0im, 0)
    end


    @testset "Commuting strings on different sites" begin
        basis = PauliStringBasis(2)

        xi = pauli_index(
            basis,
            (1, 0),
        )

        iz = pauli_index(
            basis,
            (0, 3),
        )

        @test pauli_commutator(
            basis,
            xi,
            iz,
        ) == (0im, 0)

        @test !pauli_anticommutes(
            basis,
            xi,
            iz,
        )
    end


    @testset "Anticommuting multi-site strings" begin
        basis = PauliStringBasis(2)

        xi = pauli_index(
            basis,
            (1, 0),
        )

        yz = pauli_index(
            basis,
            (2, 3),
        )

        zz = pauli_code(
            basis,
            (3, 3),
        )

        @test pauli_commutator(
            basis,
            xi,
            yz,
        ) == (2im, zz)

        @test pauli_anticommutes(
            basis,
            xi,
            yz,
        )
    end


    @testset "Two local anticommutions give global commutation" begin
        basis = PauliStringBasis(2)

        xx = pauli_index(
            basis,
            (1, 1),
        )

        yy = pauli_index(
            basis,
            (2, 2),
        )

        # XY = iZ on each site, so
        #
        # (XX)(YY) = (iZ)(iZ) = -ZZ,
        #
        # and reversing both local products gives the same global result.
        @test pauli_product(
            basis,
            xx,
            yy,
        ) == (
            -1 + 0im,
            pauli_code(basis, (3, 3)),
        )

        @test pauli_product(
            basis,
            yy,
            xx,
        ) == (
            -1 + 0im,
            pauli_code(basis, (3, 3)),
        )

        @test pauli_commutator(
            basis,
            xx,
            yy,
        ) == (0im, 0)
    end


    @testset "Product reverse consistency" begin
        for n in 1:3
            basis = PauliStringBasis(n)

            @test all(
                begin
                    phase_lr, code_lr = pauli_product(
                        basis,
                        left,
                        right,
                    )

                    phase_rl, code_rl = pauli_product(
                        basis,
                        right,
                        left,
                    )

                    code_lr == code_rl &&
                    (
                        phase_lr == phase_rl ||
                        phase_lr == -phase_rl
                    )
                end
                for left in 1:basis_size(basis)
                for right in 1:basis_size(basis)
            )
        end
    end


    @testset "Commutator coefficient structure" begin
        for n in 1:3
            basis = PauliStringBasis(n)

            @test all(
                begin
                    coefficient, code = pauli_commutator(
                        basis,
                        left,
                        right,
                    )

                    if iszero(coefficient)
                        code == 0
                    else
                        (
                            coefficient == 2im ||
                            coefficient == -2im
                        ) &&
                        1 <= code <= basis_size(basis)
                    end
                end
                for left in 1:basis_size(basis)
                for right in 1:basis_size(basis)
            )
        end
    end
end
