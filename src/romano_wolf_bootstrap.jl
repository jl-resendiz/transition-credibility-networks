"""
Romano-Wolf (2005, 2016) bootstrap in Julia.

Reads pre-computed X, y_hat, resid, inv_XtX, cluster_ids from CSV
(exported by export_bootstrap_inputs.py) and runs the Rademacher
cluster bootstrap. Compares results with Python output.

Usage:
    julia romano_wolf_bootstrap.jl
"""

using DelimitedFiles
using Random
using LinearAlgebra
using Printf
using Statistics

const B = 999
const SEED = 42
const MONTH_POSTS = [3]
const CHANNEL_VARS = ["w_geo", "w_fuel", "w_reg"]
const DATA_DIR = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "data")

# ── Load matrices for each window ──

struct WindowData
    X::Matrix{Float64}
    y_hat::Vector{Float64}
    resid::Vector{Float64}
    cluster_ids::Vector{Int}
    inv_XtX::Matrix{Float64}
    n::Int
    k::Int
    spec_vars::Vector{String}
    cluster_map::Dict{Int, Vector{Int}}
    beta_obs::Vector{Float64}
    se_obs::Vector{Float64}
end

function load_window(post::Int)::Union{WindowData, Nothing}
    x_path = joinpath(DATA_DIR, "X_$post.csv")
    isfile(x_path) || return nothing

    # X matrix: skip header row
    X_raw = readdlm(x_path, ',', Float64; header=true)[1]
    n, k = size(X_raw)

    # Vectors
    vec_raw = readdlm(joinpath(DATA_DIR, "vectors_$post.csv"), ','; header=true)[1]
    y_hat = Float64.(vec_raw[:, 2])
    resid = Float64.(vec_raw[:, 3])
    cluster_ids = Int.(vec_raw[:, 4])

    # inv_XtX
    inv_raw = readdlm(joinpath(DATA_DIR, "inv_XtX_$post.csv"), ',', Float64; header=true)[1]

    # spec_vars
    spec_vars = readlines(joinpath(DATA_DIR, "spec_vars_$post.txt"))

    # Build cluster map
    cmap = Dict{Int, Vector{Int}}()
    for i in 1:n
        cid = cluster_ids[i]
        if haskey(cmap, cid)
            push!(cmap[cid], i)
        else
            cmap[cid] = [i]
        end
    end

    # Observed betas and SEs for centred bootstrap
    obs_path = joinpath(DATA_DIR, "obs_beta_se_$post.csv")
    if isfile(obs_path)
        obs_raw = readdlm(obs_path, ',', Float64; header=true)[1]
        beta_obs = obs_raw[:, 1]
        se_obs = obs_raw[:, 2]
    else
        beta_obs = zeros(k)
        se_obs = ones(k)
    end

    return WindowData(X_raw, y_hat, resid, cluster_ids, inv_raw,
                      n, k, spec_vars, cmap, beta_obs, se_obs)
end

# ── Clustered covariance and bootstrap t-stats ──

function cluster_cov!(S::Matrix{Float64}, X::Matrix{Float64},
                      resid::Vector{Float64}, cmap::Dict{Int, Vector{Int}},
                      k::Int)
    fill!(S, 0.0)
    xu = zeros(Float64, k)
    for idxs in values(cmap)
        fill!(xu, 0.0)
        for i in idxs
            ri = resid[i]
            @inbounds for a in 1:k
                xu[a] += X[i, a] * ri
            end
        end
        @inbounds for a in 1:k
            xua = xu[a]
            for b in a:k
                v = xua * xu[b]
                S[a, b] += v
                if a != b
                    S[b, a] += v
                end
            end
        end
    end
    return S
end

function bootstrap_t_stats!(t_out::Vector{Float64}, y_star::Vector{Float64},
                            wd::WindowData,
                            S_buf::Matrix{Float64}, Xty_buf::Vector{Float64},
                            beta_buf::Vector{Float64}, resid_buf::Vector{Float64})
    X = wd.X
    inv_XtX = wd.inv_XtX
    n = wd.n
    k = wd.k

    # X'y*
    fill!(Xty_buf, 0.0)
    @inbounds for i in 1:n
        yi = y_star[i]
        for a in 1:k
            Xty_buf[a] += X[i, a] * yi
        end
    end

    # beta* = inv_XtX * X'y*
    mul!(beta_buf, inv_XtX, Xty_buf)

    # residuals
    @inbounds for i in 1:n
        fitted = 0.0
        for a in 1:k
            fitted += X[i, a] * beta_buf[a]
        end
        resid_buf[i] = y_star[i] - fitted
    end

    # Clustered covariance
    cluster_cov!(S_buf, X, resid_buf, wd.cluster_map, k)

    # V = inv_XtX * S * inv_XtX
    V = inv_XtX * S_buf * inv_XtX

    # Scale correction
    G = length(wd.cluster_map)
    if G > 1
        scale = (G / (G - 1)) * ((n - 1) / (n - k))
        V .*= scale
    end

    # Centred t-stats: (beta* - beta_obs) / se_obs
    # If beta_obs/se_obs available in wd, use them; otherwise fallback to beta*/se*
    if isdefined(wd, :beta_obs) && isdefined(wd, :se_obs)
        @inbounds for a in 1:k
            se_a = wd.se_obs[a]
            t_out[a] = se_a > 1e-15 ? (beta_buf[a] - wd.beta_obs[a]) / se_a : 0.0
        end
    else
        @inbounds for a in 1:k
            se = V[a, a] > 0 ? sqrt(V[a, a]) : 0.0
            t_out[a] = se > 1e-15 ? beta_buf[a] / se : 0.0
        end
    end
end

# ── Main ──

function main()
    println("Loading window data...")

    windows = Dict{Int, WindowData}()
    for post in MONTH_POSTS
        wd = load_window(post)
        if wd !== nothing
            windows[post] = wd
            println("  Window [-1,+$post]: N=$(wd.n), k=$(wd.k), clusters=$(length(wd.cluster_map))")
        end
    end

    # Collect all cluster IDs
    all_cids = sort(collect(union([Set(keys(wd.cluster_map)) for wd in values(windows)]...)))
    println("  Total unique clusters: $(length(all_cids))")

    # Hypothesis labels: 3 channels x 1 window = 3
    hyp_labels = Tuple{String, Int}[]
    for post in MONTH_POSTS
        for ch in CHANNEL_VARS
            push!(hyp_labels, (ch, post))
        end
    end
    n_hyp = length(hyp_labels)

    # Pre-allocate buffers per window
    bufs = Dict{Int, NamedTuple}()
    for (post, wd) in windows
        bufs[post] = (
            y_star = zeros(Float64, wd.n),
            S = zeros(Float64, wd.k, wd.k),
            Xty = zeros(Float64, wd.k),
            beta = zeros(Float64, wd.k),
            resid = zeros(Float64, wd.n),
            t_out = zeros(Float64, wd.k),
        )
    end

    # Bootstrap
    println("\nRunning $B bootstrap replications...")
    boot_t = zeros(Float64, B, n_hyp)
    boot_max_t = zeros(Float64, B)

    rng = MersenneTwister(SEED)
    t_start = time()

    for b in 1:B
        if b % 100 == 0
            elapsed = time() - t_start
            rate = b / elapsed
            @printf("  bootstrap %d/%d (%.1f iter/s)\n", b, B, rate)
        end

        # Rademacher weights per cluster
        rademacher = Dict{Int, Float64}()
        for cid in all_cids
            rademacher[cid] = rand(rng) < 0.5 ? 1.0 : -1.0
        end

        col = 0
        for post in MONTH_POSTS
            if !haskey(windows, post)
                col += length(CHANNEL_VARS)
                continue
            end

            wd = windows[post]
            bf = bufs[post]

            # y* = y_hat + w_g * resid
            @inbounds for i in 1:wd.n
                w = get(rademacher, wd.cluster_ids[i], 1.0)
                bf.y_star[i] = wd.y_hat[i] + w * wd.resid[i]
            end

            bootstrap_t_stats!(bf.t_out, bf.y_star, wd,
                               bf.S, bf.Xty, bf.beta, bf.resid)

            # Extract channel t-stats
            # Channel indices: w_geo is index 2 (after intercept),
            # w_fuel is 3, w_reg is 4 in the spec_vars order
            for ch in CHANNEL_VARS
                col += 1
                ch_idx = findfirst(==(ch), wd.spec_vars)
                if ch_idx !== nothing
                    boot_t[b, col] = abs(bf.t_out[ch_idx + 1])  # +1 for intercept
                end
            end
        end

        boot_max_t[b] = maximum(boot_t[b, :])
    end

    elapsed = time() - t_start
    @printf("\nBootstrap complete in %.2f seconds (%.1f iter/s)\n",
            elapsed, B / elapsed)

    # ── Compute p-values ──

    # We need the original t-stats. Read them from the Python output
    # or recompute from the exported data. Let's recompute.
    obs_t = zeros(Float64, n_hyp)
    col = 0
    for post in MONTH_POSTS
        if !haskey(windows, post)
            col += length(CHANNEL_VARS)
            continue
        end
        wd = windows[post]
        # OLS t-stats from the full sample
        y = wd.y_hat .+ wd.resid  # y = y_hat + resid
        Xty = wd.X' * y
        beta = wd.inv_XtX * Xty
        fitted = wd.X * beta
        res = y .- fitted
        S = zeros(Float64, wd.k, wd.k)
        cluster_cov!(S, wd.X, res, wd.cluster_map, wd.k)
        V = wd.inv_XtX * S * wd.inv_XtX
        G = length(wd.cluster_map)
        if G > 1
            scale = (G / (G - 1)) * ((wd.n - 1) / (wd.n - wd.k))
            V .*= scale
        end
        for ch in CHANNEL_VARS
            col += 1
            ch_idx = findfirst(==(ch), wd.spec_vars)
            if ch_idx !== nothing
                idx = ch_idx + 1
                se = V[idx, idx] > 0 ? sqrt(V[idx, idx]) : 0.0
                obs_t[col] = se > 1e-15 ? beta[idx] / se : 0.0
            end
        end
    end

    # Westfall-Young max-t p-values
    maxt_p = zeros(Float64, n_hyp)
    for j in 1:n_hyp
        abs_t = abs(obs_t[j])
        maxt_p[j] = count(mt -> mt >= abs_t, boot_max_t) / B
    end

    # Romano-Wolf stepdown
    order = sortperm(abs.(obs_t); rev=true)
    rw_p = zeros(Float64, n_hyp)
    remaining = collect(1:n_hyp)

    for (step, j) in enumerate(order)
        abs_t = abs(obs_t[j])
        cnt = 0
        for b in 1:B
            step_max = maximum(boot_t[b, r] for r in remaining)
            if step_max >= abs_t
                cnt += 1
            end
        end
        p_step = cnt / B
        if step > 1
            p_step = max(p_step, rw_p[order[step - 1]])
        end
        rw_p[j] = p_step
        filter!(x -> x != j, remaining)
        isempty(remaining) && break
    end

    # ── Print results ──
    println("\n" * "="^70)
    println("JULIA ROMANO-WOLF RESULTS")
    println("="^70)
    @printf("%-10s %-12s %8s %8s %8s\n",
            "Variable", "Window", "t", "Max-t p", "RW p")
    println("-"^50)
    for (j, (ch, post)) in enumerate(hyp_labels)
        @printf("%-10s [-1,+%-2d]     %8.3f %8.4f %8.4f\n",
                ch, post, obs_t[j], maxt_p[j], rw_p[j])
    end

    # ── Write comparison output ──
    out_path = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "julia_rw_results.csv")
    open(out_path, "w") do io
        println(io, "channel,window,obs_t,maxt_p,rw_p")
        for (j, (ch, post)) in enumerate(hyp_labels)
            @printf(io, "%s,%d,%.6f,%.4f,%.4f\n",
                    ch, post, obs_t[j], maxt_p[j], rw_p[j])
        end
    end
    println("\nWrote: $out_path")
end

main()
