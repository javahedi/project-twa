@testset "CTWA discrete sampling" begin

    @testset "Method selection" begin
        @test CTWA().sampling isa
              DiscreteSampling

        @test CTWA(
            sampling=GaussianSampling(),
        ).sampling isa GaussianSampling
    end


    @testset "Single-site up state" begin
        plan =
            compile_ctwa_sampling(
                Up(),
                Clustering(
                    1,
                    1,
                ),
                DiscreteSampling(),
            )

        basis = plan.basis

        x =
            pauli_index(
                basis,
                (1,),
            )

        y =
            pauli_index(
                basis,
                (2,),
            )

        z =
            pauli_index(
                basis,
                (3,),
            )

        sample =
            sample_ctwa(
                plan;
                rng=Xoshiro(1234),
            )

        @test abs(sample[x, 1]) == 1.0
        @test abs(sample[y, 1]) == 1.0
        @test sample[z, 1] == 1.0
    end


    @testset "Single-site down state" begin
        plan =
            compile_ctwa_sampling(
                Down(),
                Clustering(
                    1,
                    1,
                ),
                DiscreteSampling(),
            )

        z =
            pauli_index(
                plan.basis,
                (3,),
            )

        sample =
            sample_ctwa(
                plan;
                rng=Xoshiro(1234),
            )

        @test sample[z, 1] == -1.0
    end


    @testset "Two-site coordinates are products" begin
        plan =
            compile_ctwa_sampling(
                Up(),
                Clustering(
                    2,
                    2,
                ),
                DiscreteSampling(),
            )

        basis = plan.basis

        xi = pauli_index(basis, (1, 0))
        ix = pauli_index(basis, (0, 1))
        iy = pauli_index(basis, (0, 2))
        iz = pauli_index(basis, (0, 3))

        xx = pauli_index(basis, (1, 1))
        xy = pauli_index(basis, (1, 2))
        xz = pauli_index(basis, (1, 3))

        zi = pauli_index(basis, (3, 0))
        zz = pauli_index(basis, (3, 3))

        sample =
            sample_ctwa(
                plan;
                rng=Xoshiro(4321),
            )

        @test sample[xx, 1] ==
              sample[xi, 1] * sample[ix, 1]

        @test sample[xy, 1] ==
              sample[xi, 1] * sample[iy, 1]

        @test sample[xz, 1] ==
              sample[xi, 1] * sample[iz, 1]

        @test sample[zi, 1] == 1.0
        @test sample[iz, 1] == 1.0
        @test sample[zz, 1] == 1.0
    end


    @testset "Float32" begin
        plan =
            compile_ctwa_sampling(
                Up(),
                Clustering(
                    2,
                    2,
                ),
                DiscreteSampling();
                T=Float32,
            )

        sample =
            sample_ctwa(
                plan;
                rng=Xoshiro(99),
            )

        @test eltype(sample) == Float32
    end
end