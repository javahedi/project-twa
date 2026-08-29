@testset "CTWA sampling" begin

    @testset "Single-site up state" begin
        plan =
            compile_ctwa_sampling(
                Up(),
                Clustering(
                    1,
                    1,
                ),
            )

        basis = plan.basis

        xgen =
            pauli_index(
                basis,
                (1,),
            )

        ygen =
            pauli_index(
                basis,
                (2,),
            )

        zgen =
            pauli_index(
                basis,
                (3,),
            )

        @test plan.means[xgen, 1] == 0.0
        @test plan.means[ygen, 1] == 0.0
        @test plan.means[zgen, 1] == 1.0

        @test plan.flip_masks[xgen, 1] != 0
        @test plan.flip_masks[ygen, 1] != 0
        @test plan.flip_masks[zgen, 1] == 0
    end


    @testset "Single-site down state" begin
        plan =
            compile_ctwa_sampling(
                Down(),
                Clustering(
                    1,
                    1,
                ),
            )

        basis = plan.basis

        zgen =
            pauli_index(
                basis,
                (3,),
            )

        @test plan.means[zgen, 1] == -1.0
    end


    @testset "Two-site up cluster exact means" begin
        plan =
            compile_ctwa_sampling(
                Up(),
                Clustering(
                    2,
                    2,
                ),
            )

        basis = plan.basis

        zi =
            pauli_index(
                basis,
                (3, 0),
            )

        iz =
            pauli_index(
                basis,
                (0, 3),
            )

        zz =
            pauli_index(
                basis,
                (3, 3),
            )

        xx =
            pauli_index(
                basis,
                (1, 1),
            )

        @test plan.means[zi, 1] == 1.0
        @test plan.means[iz, 1] == 1.0
        @test plan.means[zz, 1] == 1.0
        @test plan.means[xx, 1] == 0.0
    end


    @testset "Shared flip masks reproduce exact correlations" begin
        plan =
            compile_ctwa_sampling(
                Up(),
                Clustering(
                    2,
                    2,
                ),
            )

        basis = plan.basis

        ix =
            pauli_index(
                basis,
                (0, 1),
            )

        iy =
            pauli_index(
                basis,
                (0, 2),
            )

        zx =
            pauli_index(
                basis,
                (3, 1),
            )

        zy =
            pauli_index(
                basis,
                (3, 2),
            )

        xx =
            pauli_index(
                basis,
                (1, 1),
            )

        yy =
            pauli_index(
                basis,
                (2, 2),
            )

        # IX and ZX act on |↑↑> with the same flipped product state and
        # identical phase, so their sampled fluctuations must be identical.
        @test plan.flip_masks[ix, 1] ==
              plan.flip_masks[zx, 1]

        @test plan.amplitude_re[ix, 1] ==
              plan.amplitude_re[zx, 1]

        @test plan.amplitude_im[ix, 1] ==
              plan.amplitude_im[zx, 1]

        # Likewise for IY and ZY.
        @test plan.flip_masks[iy, 1] ==
              plan.flip_masks[zy, 1]

        # XX and YY flip both sites, but YY carries i*i = -1 relative to XX.
        @test plan.flip_masks[xx, 1] ==
              plan.flip_masks[yy, 1]

        @test plan.amplitude_re[xx, 1] == 1.0
        @test plan.amplitude_re[yy, 1] == -1.0
    end


    @testset "Sample obeys deterministic cluster relations" begin
        plan =
            compile_ctwa_sampling(
                Up(),
                Clustering(
                    2,
                    2,
                ),
            )

        basis = plan.basis

        sample =
            sample_ctwa(
                plan;
                rng=Xoshiro(1234),
            )

        ix =
            pauli_index(
                basis,
                (0, 1),
            )

        zx =
            pauli_index(
                basis,
                (3, 1),
            )

        iy =
            pauli_index(
                basis,
                (0, 2),
            )

        zy =
            pauli_index(
                basis,
                (3, 2),
            )

        xx =
            pauli_index(
                basis,
                (1, 1),
            )

        yy =
            pauli_index(
                basis,
                (2, 2),
            )

        zz =
            pauli_index(
                basis,
                (3, 3),
            )

        @test sample[ix, 1] == sample[zx, 1]
        @test sample[iy, 1] == sample[zy, 1]
        @test sample[yy, 1] == -sample[xx, 1]
        @test sample[zz, 1] == 1.0
    end


    @testset "Domain wall cluster means" begin
        state =
            DomainWall()

        clustering =
            Clustering(
                4,
                2,
            )

        plan =
            compile_ctwa_sampling(
                state,
                clustering,
            )

        basis = plan.basis

        zz =
            pauli_index(
                basis,
                (3, 3),
            )

        zi =
            pauli_index(
                basis,
                (3, 0),
            )

        # First cluster is ↑↑ and second cluster is ↓↓ for N=4.
        @test plan.means[zz, 1] == 1.0
        @test plan.means[zz, 2] == 1.0

        @test plan.means[zi, 1] == 1.0
        @test plan.means[zi, 2] == -1.0
    end


    @testset "Gaussian first and second moments for one spin" begin
        plan =
            compile_ctwa_sampling(
                Up(),
                Clustering(
                    1,
                    1,
                ),
            )

        rng =
            Xoshiro(4321)

        nsamples = 20_000

        sx = 0.0
        sy = 0.0
        sx2 = 0.0
        sy2 = 0.0
        sxy = 0.0

        for _ in 1:nsamples
            sample =
                sample_ctwa(
                    plan;
                    rng=rng,
                )

            x = sample[1, 1]
            y = sample[2, 1]

            sx += x
            sy += y
            sx2 += x * x
            sy2 += y * y
            sxy += x * y
        end

        mx = sx / nsamples
        my = sy / nsamples

        vx =
            sx2 / nsamples -
            mx^2

        vy =
            sy2 / nsamples -
            my^2

        cxy =
            sxy / nsamples -
            mx * my

        @test abs(mx) < 0.03
        @test abs(my) < 0.03

        @test vx ≈ 1.0 atol=0.04
        @test vy ≈ 1.0 atol=0.04
        @test abs(cxy) < 0.04
    end


    @testset "Float32 sampling" begin
        plan =
            compile_ctwa_sampling(
                Up(),
                Clustering(
                    2,
                    1,
                );
                T=Float32,
            )

        sample =
            sample_ctwa(
                plan;
                rng=Xoshiro(99),
            )

        @test eltype(sample) == Float32
    end


    @testset "Preallocated sampling" begin
        plan =
            compile_ctwa_sampling(
                Up(),
                Clustering(
                    2,
                    1,
                ),
            )

        x =
            zeros(
                Float64,
                size(plan.means),
            )

        result =
            sample_ctwa!(
                x,
                plan;
                rng=Xoshiro(7),
            )

        @test result === x

        wrong =
            zeros(
                Float64,
                2,
                2,
            )

        @test_throws DimensionMismatch sample_ctwa!(
            wrong,
            plan;
            rng=Xoshiro(7),
        )
    end
end
