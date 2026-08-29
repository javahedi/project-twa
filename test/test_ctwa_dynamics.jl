@testset "CTWA deterministic dynamics" begin

    @testset "Single spin in z field" begin
        model = SpinModel(
            Chain(1),
            Field(
                :z,
                0.7,
            ),
        )

        clustering = Clustering(
            1,
            1,
        )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )

        # Basis ordering for one spin is X,Y,Z.
        x = reshape(
            [
                0.4,
                -0.2,
                0.9,
            ],
            3,
            1,
        )

        dx = similar(x)

        ctwa_rhs!(
            dx,
            x,
            compiled,
            0.0,
        )

        # Pauli convention:
        #
        #     ds/dt = 2 h × s.
        #
        # For h=(0,0,hz):
        #
        #     dx/dt = -2 hz y
        #     dy/dt =  2 hz x
        #     dz/dt =  0.
        @test dx[1, 1] ≈ -2 * 0.7 * x[2, 1]
        @test dx[2, 1] ≈  2 * 0.7 * x[1, 1]
        @test dx[3, 1] ≈ 0.0
    end


    @testset "Single-site clusters reproduce DTWA RHS" begin
        model = SpinModel(
            Chain(
                4,
                PeriodicBoundary(),
            ),
            XXZ(
                J=0.8,
                Δ=0.6,
            ) +
            Field(
                :x,
                0.2,
            ) +
            Field(
                :z,
                -0.35,
            ),
        )

        clustering = Clustering(
            4,
            1,
        )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )

        # Same physical phase-space point in the two layouts:
        #
        # DTWA: (site, component)
        # CTWA: (generator, cluster)
        u_dtwa = [
             0.2   0.7  -0.4
            -0.8   0.1   0.5
             0.6  -0.3   0.9
            -0.1  -0.9   0.2
        ]

        x_ctwa = permutedims(
            u_dtwa,
            (2, 1),
        )

        du_dtwa = similar(
            u_dtwa,
        )

        dx_ctwa = similar(
            x_ctwa,
        )

        dtwa_rhs!(
            du_dtwa,
            u_dtwa,
            model,
            0.0,
        )

        ctwa_rhs!(
            dx_ctwa,
            x_ctwa,
            compiled,
            0.0,
        )

        @test dx_ctwa ≈ permutedims(
            du_dtwa,
            (2, 1),
        )
    end


    @testset "Two-site cluster internal interaction" begin
        model = SpinModel(
            Chain(2),
            XXZ(
                J=1.0,
                Δ=0.5,
            ),
        )

        clustering = Clustering(
            2,
            2,
        )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )

        basis = compiled.basis

        x = zeros(
            Float64,
            basis_size(basis),
            1,
        )

        # Give the phase-space point several nonzero coordinates so the
        # intra-cluster XX, YY, ZZ Hamiltonian has something to rotate.
        x[
            pauli_index(basis, (1, 0)),
            1,
        ] = 0.3

        x[
            pauli_index(basis, (2, 0)),
            1,
        ] = -0.4

        x[
            pauli_index(basis, (0, 3)),
            1,
        ] = 0.8

        x[
            pauli_index(basis, (1, 2)),
            1,
        ] = -0.2

        dx = similar(x)

        ctwa_rhs!(
            dx,
            x,
            compiled,
            0.0,
        )

        @test all(
            isfinite,
            dx,
        )

        @test any(
            !iszero,
            dx,
        )
    end


    @testset "Inter-cluster interaction affects both clusters" begin
        model = SpinModel(
            Chain(4),
            XXZ(
                J=1.0,
                Δ=1.0,
            ),
        )

        clustering = Clustering(
            4,
            2,
        )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )

        basis = compiled.basis

        x = zeros(
            Float64,
            basis_size(basis),
            2,
        )

        # Populate a generic deterministic phase-space point.
        for cluster in 1:2
            for generator in 1:basis_size(basis)
                x[
                    generator,
                    cluster,
                ] =
                    0.01 *
                    (
                        generator +
                        3 * cluster
                    )
            end
        end

        dx = similar(x)

        ctwa_rhs!(
            dx,
            x,
            compiled,
            0.0,
        )

        @test any(
            !iszero,
            view(dx, :, 1),
        )

        @test any(
            !iszero,
            view(dx, :, 2),
        )
    end


    @testset "Zero Hamiltonian gives zero RHS" begin
        model = SpinModel(
            Chain(2),
            Field(
                :z,
                0.0,
            ),
        )

        clustering = Clustering(
            2,
            1,
        )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )

        x = [
             0.2   0.5
            -0.1   0.3
             0.7  -0.8
        ]

        dx = similar(x)

        ctwa_rhs!(
            dx,
            x,
            compiled,
            0.0,
        )

        @test iszero(
            dx,
        )
    end


    @testset "Float32 propagation" begin
        model = SpinModel(
            Chain(2),
            XXZ(
                J=Float32(1),
                Δ=Float32(0.5),
            ),
        )

        clustering = Clustering(
            2,
            1,
        )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering;
                T=Float32,
            )

        x = Float32[
             0.2   0.4
            -0.3   0.5
             0.8  -0.1
        ]

        dx = similar(x)

        ctwa_rhs!(
            dx,
            x,
            compiled,
            Float32(0),
        )

        @test eltype(dx) == Float32
        @test all(
            isfinite,
            dx,
        )
    end


    @testset "State-shape validation" begin
        model = SpinModel(
            Chain(4),
            XXZ(),
        )

        clustering = Clustering(
            4,
            2,
        )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )

        dimension =
            basis_size(
                compiled.basis,
            )

        correct = zeros(
            dimension,
            2,
        )

        wrong_generators = zeros(
            dimension - 1,
            2,
        )

        wrong_clusters = zeros(
            dimension,
            1,
        )

        @test_throws DimensionMismatch ctwa_rhs!(
            wrong_generators,
            correct,
            compiled,
            0.0,
        )

        @test_throws DimensionMismatch ctwa_rhs!(
            correct,
            wrong_clusters,
            compiled,
            0.0,
        )
    end
end
