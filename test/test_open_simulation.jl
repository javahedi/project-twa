using Random

@testset "Open DTWA and CTWA" begin
    @test_throws ArgumentError LocalDissipation(decay=-1)
    @test_throws ArgumentError LocalDissipation(dephasing=Inf)

    model = SpinModel(Chain(1), XXZ(J=0.0))
    dtwa = DTWA(trajectories=8)
    ctwa = CTWA(cluster_size=1, trajectories=8)
    rates = LocalDissipation(decay=0.4, dephasing=0.3)
    for method in (dtwa, ctwa)
        closed = simulate(model, Up(), method; tspan=(0.0, 0.5), rng=Xoshiro(1))
        zero_rates = simulate(
            model, Up(), method;
            tspan=(0.0, 0.5), rng=Xoshiro(1), dissipation=LocalDissipation(),
        )
        @test closed.trajectories == zero_rates.trajectories
        @test_throws ArgumentError simulate(
            model, Up(), method; tspan=(0.0, 0.5), dissipation=rates, dt=0,
        )
    end

    options = (; tspan=(0.0, 0.5), saveat=[0.0, 0.23, 0.5], dt=0.01, dissipation=rates)
    single = simulate(model, Up(), dtwa; options..., rng=Xoshiro(12))
    cluster = simulate(model, Up(), ctwa; options..., rng=Xoshiro(12))
    @test times(single) == times(cluster) == options.saveat
    @test single.trajectories == permutedims(cluster.trajectories, (2, 1, 3, 4))

    for method in (DTWA(trajectories=1200), CTWA(cluster_size=1, trajectories=1200))
        decay = simulate(
            model, Up(), method;
            tspan=(0.0, 0.5), saveat=[0.5], dt=0.01, rng=Xoshiro(2),
            dissipation=LocalDissipation(decay=0.4),
        )
        @test only(expectation(decay, LocalMagnetization(:z, 1))) ≈ 2exp(-0.2) - 1 atol=0.05

        dephasing = simulate(
            model, Polarized(:x), method;
            tspan=(0.0, 0.5), saveat=[0.5], dt=0.01, rng=Xoshiro(3),
            dissipation=LocalDissipation(dephasing=0.6),
        )
        @test only(expectation(dephasing, LocalMagnetization(:x, 1))) ≈ exp(-0.3) atol=0.05
    end

    pair = simulate(
        SpinModel(Chain(2), XXZ(J=0.0)), Up(), CTWA(trajectories=1500);
        tspan=(0.0, 0.5), saveat=[0.5], dt=0.01, rng=Xoshiro(4),
        dissipation=rates,
    )
    @test only(expectation(pair, Correlation(:z, 1, :z, 2))) ≈ (2exp(-0.2) - 1)^2 atol=0.06

    # In the one-excitation sector, XX exchange oscillates at 4J.
    # Emission transfers the excitation to the all-down state.
    interacting = simulate(
        SpinModel(Chain(2), XXZ(J=0.3, Δ=0.0)),
        DomainWall(),
        CTWA(trajectories=3000);
        tspan=(0.0, 0.5), saveat=[0.5], dt=0.005, rng=Xoshiro(5),
        dissipation=LocalDissipation(decay=0.4),
    )
    z1 = -1 + exp(-0.2) * (1 + cos(0.6))
    @test only(expectation(interacting, LocalMagnetization(:z, 1))) ≈ z1 atol=0.06

end
