@testset "Compiled CTWA algebra" begin

    @testset "One-spin action contents" begin
        model = SpinModel(
            Chain(1),
            Field(
                :z,
                1.0,
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

        algebra =
            compile_ctwa_algebra(
                compiled,
            )

        basis = compiled.basis

        x = pauli_index(
            basis,
            (1,),
        )

        y = pauli_index(
            basis,
            (2,),
        )

        z = pauli_index(
            basis,
            (3,),
        )

        action =
            commutator_action(
                algebra,
                z,
            )

        entries = Dict(
            entry.state_generator => (
                entry.result_generator,
                entry.sign,
            )
            for entry in action.entries
        )

        @test length(entries) == 2

        # [X,Z] = -2i Y
        @test entries[x] == (
            y,
            Int8(-1),
        )

        # [Y,Z] = +2i X
        @test entries[y] == (
            x,
            Int8(1),
        )
    end


    @testset "Only used generators are cached" begin
        model = SpinModel(
            Chain(2),
            Field(
                :z,
                0.5,
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

        algebra =
            compile_ctwa_algebra(
                compiled,
            )

        basis = compiled.basis

        x = pauli_index(
            basis,
            (1,),
        )

        z = pauli_index(
            basis,
            (3,),
        )

        @test algebra.actions[z] !== nothing
        @test algebra.actions[x] === nothing

        @test_throws ArgumentError commutator_action(
            algebra,
            x,
        )
    end


    @testset "Cached RHS matches reference RHS" begin
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

        for cluster_size in (
            1,
            2,
            4,
        )
            clustering =
                Clustering(
                    4,
                    cluster_size,
                )

            compiled =
                compile_ctwa_hamiltonian(
                    model,
                    clustering,
                )

            algebra =
                compile_ctwa_algebra(
                    compiled,
                )

            dimension =
                basis_size(
                    compiled.basis,
                )

            nclusters =
                cluster_count(
                    clustering,
                )

            x = reshape(
                [
                    0.01 * (
                        generator +
                        7 * cluster
                    )
                    for generator in 1:dimension
                    for cluster in 1:nclusters
                ],
                dimension,
                nclusters,
            )

            dx_reference =
                similar(x)

            dx_cached =
                similar(x)

            ctwa_rhs!(
                dx_reference,
                x,
                compiled,
                0.0,
            )

            ctwa_rhs!(
                dx_cached,
                x,
                compiled,
                algebra,
                0.0,
            )

            @test dx_cached ≈ dx_reference
        end
    end


    @testset "Cluster-size-one cached RHS matches DTWA" begin
        model = SpinModel(
            Chain(4),
            XXZ(
                J=1.1,
                Δ=0.7,
            ) +
            Field(
                :y,
                0.15,
            ),
        )

        clustering =
            Clustering(
                4,
                1,
            )

        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )

        algebra =
            compile_ctwa_algebra(
                compiled,
            )

        u_dtwa = [
             0.2   0.7  -0.4
            -0.8   0.1   0.5
             0.6  -0.3   0.9
            -0.1  -0.9   0.2
        ]

        x_ctwa =
            permutedims(
                u_dtwa,
                (2, 1),
            )

        du_dtwa =
            similar(
                u_dtwa,
            )

        dx_ctwa =
            similar(
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
            algebra,
            0.0,
        )

        @test dx_ctwa ≈
            permutedims(
                du_dtwa,
                (2, 1),
            )
    end
end
