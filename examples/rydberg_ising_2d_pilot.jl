# 2D Rydberg-Ising pilot: dTWA vs spatial dcTWA
#
# Inspired by the 2D transverse-field Ising dynamics studied by
# Bensch et al., arXiv:2608.07211 (2026).
#
# This is a small package-validation example, not a reproduction of the
# experimental geometry or full paper analysis.
#
# Package convention:
#
#   H = sum_{i<j,a} J^a_ij sigma_i^a sigma_j^a
#       + sum_{i,a} h_i^a sigma_i^a
#
# The Rydberg-Ising Hamiltonian
#
#   H / ħ = (1/8) sum_{i != j} U_ij Z_i Z_j
#           + (Ω/2) sum_i X_i
#
# therefore maps to J_ij^z = U_ij / 4 when the package sums over i < j.
using Plots
using TWA
using Random
using Printf

const U0 = 1.0

# Keep the pilot deliberately small.
const LATTICE = SquareLattice(8, 8, PeriodicBoundary())
const STATE = Down()

const OMEGAS = (0.8, 1.2, 1.6)
const TMAX = 3.0
const DT = 0.05
const NTRAJ = 32
const BASE_SEED = 20260905

"""
    rydberg_ising_model(Ω)

Clean periodic 2D transverse-field Ising model with van der Waals
power-law interactions U(r) = U0 / r^6 and h_z = 0.
"""
function rydberg_ising_model(Ω)
    return SpinModel(
        LATTICE,
        PowerLaw(
            Ising(:z; J=U0 / 4);
            α=6.0,
        ) +
        Field(:x, Ω / 2),
    )
end

"""
    methods_for_pilot()

Return the three discrete phase-space approximations used in this pilot:

- dTWA: independent single-spin trajectories;
- dcTWA k=2: spatial 2×1 clusters;
- dcTWA k=4: spatial 2×2 clusters.
"""
function methods_for_pilot()
    c2 = Clustering(LATTICE, (2, 1))  # k = 2
    c4 = Clustering(LATTICE, (2, 2))  # k = 4
    c8 = Clustering(LATTICE, (4, 2))  # k = 8

    return (
        ("dTWA", DTWA(trajectories=NTRAJ)),
        (
            "dcTWA k=2",
            CTWA(
                clustering=c2,
                trajectories=NTRAJ,
                sampling=DiscreteSampling(),
            ),
        ),
        (
            "dcTWA k=4",
            CTWA(
                clustering=c4,
                trajectories=NTRAJ,
                sampling=DiscreteSampling(),
            ),
        ),

        

        (
            "dcTWA k=8",
            CTWA(
                clustering=c8,
                trajectories=NTRAJ,
                sampling=DiscreteSampling(),
            ),
        ),
    )
end

"""
    run_pilot()

Run all Ω values and return a nested dictionary

    results[Ω][label] = (result=result, mz=mz)

where `mz` is the magnetization density ⟨S_z⟩ / N.
"""
function run_pilot()
    results = Dict{Float64, Dict{String, NamedTuple}}()

    for (ω_index, Ω) in enumerate(OMEGAS)
        @printf("\nΩ/U0 = %.2f\n", Ω / U0)

        model = rydberg_ising_model(Ω)
        per_method = Dict{String, NamedTuple}()

        for (method_index, (label, method)) in enumerate(methods_for_pilot())
            seed = BASE_SEED + 100 * ω_index + method_index
            rng = MersenneTwister(seed)

            @printf("  %-12s ... ", label)
            elapsed = @elapsed begin
                result = simulate(
                    model,
                    STATE,
                    method;
                    tspan=(0.0, TMAX),
                    saveat=DT,
                    rng=rng,
                )

                mz = expectation(
                    result,
                    Magnetization(:z),
                )

                per_method[label] = (
                    result=result,
                    mz=mz,
                )
            end

            @printf(
                "done in %.2f s, Mz(tmax) = % .6f\n",
                elapsed,
                last(per_method[label].mz),
            )
        end

        results[Float64(Ω)] = per_method
    end

    return results
end

"""
    write_csv(results, filename="rydberg_ising_2d_pilot.csv")

Write the magnetization curves to a simple CSV file with columns

    omega, method, time, mz
"""
function write_csv(
    results,
    filename="rydberg_ising_2d_pilot.csv",
)
    open(filename, "w") do io
        println(io, "omega,method,time,mz")

        for Ω in sort!(collect(keys(results)))
            for (label, data) in results[Ω]
                t = times(data.result)

                for index in eachindex(t, data.mz)
                    println(
                        io,
                        Ω, ",",
                        label, ",",
                        t[index], ",",
                        data.mz[index],
                    )
                end
            end
        end
    end

    return filename
end


function plot_pilot(
    results;
    filename="rydberg_ising_2d_pilot.png",
)
    panels = []

    for Ω in OMEGAS
        p = plot(
            xlabel="U₀ t",
            ylabel="⟨Sᶻ⟩ / N",
            title="Ω/U₀ = $(Ω / U0)",
            legend=:bottomright,
            linewidth=2,
            framestyle=:box,
        )

        for label in ("dTWA", "dcTWA k=2", "dcTWA k=4")
            data = results[Float64(Ω)]

            plot!(
                p,
                times(data[label].result),
                data[label].mz;
                label=label,
            )
        end

        push!(panels, p)
    end

    figure = plot(
        panels...;
        layout=(1, 3),
        size=(1200, 360),
        margin=4Plots.mm,
    )

    savefig(figure, filename)

    return figure
end


if abspath(PROGRAM_FILE) == @__FILE__
    println("2D Rydberg-Ising pilot")
    println("=======================")
    println("lattice      = 8 × 8 periodic square")
    println("U0           = ", U0)
    println("Ω/U0         = ", OMEGAS)
    println("trajectories = ", NTRAJ)
    println("tmax         = ", TMAX)
    println("saveat       = ", DT)

    results = run_pilot()

    csv_filename = write_csv(results)
    plot_pilot(results)

    println("\nSaved magnetization curves to: ", csv_filename)
    println("Saved figure to: rydberg_ising_2d_pilot.png")


    
end
