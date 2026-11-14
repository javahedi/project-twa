# Reproduce the parameter sets in legacy example_staggered_longrange_open.jl.
# julia --project examples/open_twa_legacy_neel.jl 2.0 0.0 decay 2 [disorders] [trajectories] [output.csv] [N6|N12]
using Random
using TWA

struct DilutedChain <: TWA.AbstractGeometry
    positions::Vector{Int}
end

TWA.nsites(chain::DilutedChain) = length(chain.positions)
TWA.distance(chain::DilutedChain, i::Integer, j::Integer) =
    abs(chain.positions[i] - chain.positions[j])
function TWA._check_site(chain::DilutedChain, site::Integer)
    checkbounds(chain.positions, site)
    return nothing
end

struct LatticeNeel <: AbstractState
    signs::Vector{Int}
end

TWA.state_direction(state::LatticeNeel, site::Integer, ::Integer) =
    (0, 0, state.signs[site])

function occupied_positions(nsites, lattice_length, rng)
    # These legacy presets use StatsBase.sample!'s partial Fisher–Yates path.
    candidates = collect(1:lattice_length)
    for i in 1:nsites
        j = rand(rng, i:lattice_length)
        candidates[i], candidates[j] = candidates[j], candidates[i]
    end
    return sort!(candidates[1:nsites])
end

function strongest_bond_clustering(model, cluster_size)
    n = nsites(geometry(model))
    n % cluster_size == 0 ||
        throw(ArgumentError("cluster size must divide the number of occupied sites"))

    strength(i, j) = maximum(abs, coupling_components(model, i, j))
    remaining = collect(1:n)
    clusters = Vector{Int}[]
    while !isempty(remaining)
        pairs = [(i, j) for i in remaining for j in remaining if i < j]
        best_pair = argmax(pair -> strength(pair...), pairs)
        cluster = collect(best_pair)
        while length(cluster) < cluster_size
            candidates = setdiff(remaining, cluster)
            site = argmax(i -> sum(j -> strength(i, j), cluster), candidates)
            push!(cluster, site)
        end
        push!(clusters, cluster)
        remaining = setdiff(remaining, cluster)
    end
    return Clustering(n, reduce(hcat, clusters; init=Matrix{Int}(undef, cluster_size, 0)))
end

function staggered_magnetization(result, signs)
    return sum(sign * expectation(result, LocalMagnetization(:z, i))
               for (i, sign) in enumerate(signs)) / length(signs)
end

function main(args)
    4 <= length(args) <= 8 ||
        error("Usage: ALPHA DELTA decay|dephasing CLUSTER_SIZE [DISORDERS] [TRAJECTORIES] [OUTPUT.csv] [N6|N12]")
    α, Δ = parse.(Float64, args[1:2])
    channel = args[3]
    channel in ("decay", "dephasing") || error("channel must be decay or dephasing")
    cluster_size = parse(Int, args[4])
    preset = get(args, 8, "N6")
    n, lattice_length, rate, default_disorders = if preset == "N6"
        (6, 60, 0.1, 200)
    elseif preset == "N12"
        (12, 120, 0.5, 500)
    else
        error("preset must be N6 or N12")
    end
    cluster_size in (2, 3, 4) && n % cluster_size == 0 ||
        error("choose a cluster size of 2, 3, or 4 that divides N")
    disorders = length(args) >= 5 ? parse(Int, args[5]) : default_disorders
    trajectories = length(args) >= 6 ? parse(Int, args[6]) : 100
    disorders > 0 || error("disorders must be positive")
    output = get(args, 7, "open_twa_legacy_neel.csv")
    saveat = 10 .^ range(-1.0, 2.0; length=80)
    dissipation = channel == "decay" ? LocalDissipation(decay=rate) :
                                      LocalDissipation(dephasing=rate)
    dtwa = zeros(length(saveat))
    ctwa = zeros(length(saveat))
    dtwa_squared = zeros(length(saveat))
    ctwa_squared = zeros(length(saveat))

    dtwa_method = DTWA(; trajectories)
    for disorder in 1:disorders
        positions = occupied_positions(n, lattice_length, MersenneTwister(1234 + disorder))
        signs = [isodd(position) ? 1 : -1 for position in positions]
        state = LatticeNeel(signs)
        model = SpinModel(DilutedChain(positions), PowerLaw(XXZ(J=1.0, Δ=Δ); α))
        clustering = strongest_bond_clustering(model, cluster_size)
        ctwa_method = CTWA(; clustering, trajectories)
        for (method, total, squared) in (
            (dtwa_method, dtwa, dtwa_squared), (ctwa_method, ctwa, ctwa_squared),
        )
            result = simulate(
                model, state, method;
                dissipation, saveat, dt=5e-3,
                tspan=(first(saveat), last(saveat)),
                rng=Xoshiro(1234 + 10_000 + disorder),
            )
            values = staggered_magnetization(result, signs)
            total .+= values
            squared .+= abs2.(values)
        end
    end

    open(output, "w") do io
        println(io, "time,open_dtwa_mstz,open_dtwa_sem,open_ctwa_mstz,open_ctwa_sem")
        for i in eachindex(saveat)
            dtwa_mean = dtwa[i] / disorders
            ctwa_mean = ctwa[i] / disorders
            dtwa_sem = disorder_sem(dtwa_mean, dtwa_squared[i], disorders)
            ctwa_sem = disorder_sem(ctwa_mean, ctwa_squared[i], disorders)
            println(io, join((saveat[i], dtwa_mean, dtwa_sem, ctwa_mean, ctwa_sem), ','))
        end
    end
    println("Wrote ", output)
end

function disorder_sem(mean, squared, count)
    return count == 1 ? NaN : sqrt(max(0, squared / count - mean^2) / (count - 1))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
