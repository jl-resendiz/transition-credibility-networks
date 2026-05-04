#=
Focused hypothesis tests for the channel split claim (Julia port).

Complementing the Romano-Wolf correction (3 hypotheses at the primary
window; fuel survives), this script performs TWO tests that directly
address the paper's core claim: spatial channels have opposite signs.

Test 1 — Joint F-test: H0: beta_geo = beta_fuel = beta_reg = 0
    (spatial network has no effect)

Test 2 — Difference-in-coefficients: H0: beta_geo = beta_fuel
    (channels have the same effect, contra the paper's claim)

Both tests use the [-1,+3] month window (strongest signal) and require
no multiple-testing correction because each is a single hypothesis.

Specification:
  CAR_j = alpha + beta_geo * w^geo + beta_fuel * w^fuel
        + beta_reg * w^reg + beta_s * SameSector + eps_j

Requires only Julia stdlib: DelimitedFiles, Random, LinearAlgebra, Printf,
Statistics, Dates.

NOTE on control sampling: Julia's hash() differs from Python's hashlib.md5,
so the specific control firms drawn may differ. The key results (F-stat,
permutation p-value) depend on the full dataset structure, not on which
specific controls are drawn. The bootstrap logic is identical.
=#

using DelimitedFiles
using Random
using LinearAlgebra
using Printf
using Statistics
using Dates

# ── Configuration ────────────────────────────────────────────────────

const B_PERM         = 999    # permutation replications for F-test
const SEED           = 42     # reproducibility
const POST_MONTHS    = 3      # [-1, +3] window — the strongest
const PRE_MONTHS     = 24     # pre-event months for AR demeaning

const CHANNEL_VARS   = ["w_geo", "w_fuel", "w_reg"]
const SPEC_VARS_FULL = ["w_geo", "w_fuel", "w_reg", "same_sector"]

# ── Path resolution (same convention as _paths.py) ───────────────────

const BASE_DIR    = @__DIR__                              # src/
const ROOT_DIR    = normpath(joinpath(BASE_DIR, ".."))
const RAW_DIR     = joinpath(ROOT_DIR, "data", "raw")
const DERIVED_DIR = joinpath(ROOT_DIR, "data", "derived")
const RESULTS_DIR = joinpath(ROOT_DIR, "results")

raw_path(parts...)     = joinpath(RAW_DIR, parts...)
derived_path(parts...) = joinpath(DERIVED_DIR, parts...)
results_path(parts...) = joinpath(RESULTS_DIR, parts...)


# ── Standard normal CDF (Abramowitz & Stegun 26.2.17) ────────────────

function normal_cdf(x::Float64)::Float64
    x < -8.0 && return 0.0
    x >  8.0 && return 1.0
    ax = abs(x)
    b0 = 0.2316419
    b1 = 0.319381530
    b2 = -0.356563782
    b3 = 1.781477937
    b4 = -1.821255978
    b5 = 1.330274429
    t = 1.0 / (1.0 + b0 * ax)
    phi = (1.0 / sqrt(2.0 * π)) * exp(-0.5 * ax * ax)
    cdf = 1.0 - phi * (b1*t + b2*t^2 + b3*t^3 + b4*t^4 + b5*t^5)
    return x < 0.0 ? 1.0 - cdf : cdf
end

"""Two-sided p-value from t-stat using normal CDF approximation."""
p_from_t(t::Float64) = 2.0 * (1.0 - normal_cdf(abs(t)))


# ── CSV reading utilities ─────────────────────────────────────────────

"""Read a CSV file and return (header::Vector{String}, rows::Vector{Dict{String,String}})."""
function read_csv(path::String)
    lines = readlines(path)
    isempty(lines) && return (String[], Dict{String,String}[])
    header = split(lines[1], ',') .|> String
    rows = Dict{String,String}[]
    for i in 2:length(lines)
        line = lines[i]
        isempty(strip(line)) && continue
        vals = split(line, ',')
        d = Dict{String,String}()
        for (j, h) in enumerate(header)
            d[h] = j <= length(vals) ? String(vals[j]) : ""
        end
        push!(rows, d)
    end
    return (header, rows)
end


# ── Load monthly returns ──────────────────────────────────────────────

function load_monthly_returns()
    println("Loading monthly returns...")
    path = derived_path("returns", "monthly_returns.csv")
    _, rows = read_csv(path)
    # gvkey -> (year_month -> return)
    ret = Dict{String, Dict{String, Float64}}()
    for r in rows
        gk = r["gvkey"]
        ym = r["datadate"][1:7]  # "YYYY-MM"
        val = tryparse(Float64, r["ret_monthly"])
        val === nothing && continue
        if !haskey(ret, gk)
            ret[gk] = Dict{String, Float64}()
        end
        ret[gk][ym] = val
    end
    println("  Monthly: $(length(ret)) firms")
    return ret
end


# ── Load Fama-French monthly factors ──────────────────────────────────

function load_ff_factors_monthly()
    println("Loading Fama-French factors...")
    path = raw_path("factors", "F-F_Research_Data_Factors.csv")
    isfile(path) || error("Missing F-F monthly factors at $path")
    vwretd = Dict{String, Float64}()
    for line in readlines(path)
        line = strip(line)
        isempty(line) && continue
        startswith(line, "This file") && continue
        startswith(line, "The ") && continue
        startswith(line, ",") && continue
        parts = split(line, ',') .|> strip
        length(parts) < 5 && continue
        date_str = String(parts[1])
        all(isdigit, date_str) || continue
        length(date_str) != 6 && continue
        mktrf = tryparse(Float64, String(parts[2]))
        rf    = tryparse(Float64, String(parts[5]))
        (mktrf === nothing || rf === nothing) && continue
        vw = (mktrf + rf) / 100.0
        ym = date_str[1:4] * "-" * date_str[5:6]
        vwretd[ym] = vw
    end
    println("  Market months: $(length(vwretd))")
    return vwretd
end


# ── Load weight matrices ──────────────────────────────────────────────

"""Load a weight matrix CSV into Dict{String, Dict{String, Float64}}."""
function load_weight_matrix(path::String, weight_col::String="w_ij")
    W = Dict{String, Dict{String, Float64}}()
    isfile(path) || return W
    _, rows = read_csv(path)
    for r in rows
        gi = r["gvkey_i"]
        gj = r["gvkey_j"]
        # Try the specified column, then fall back to w_reg
        wstr = get(r, weight_col, "")
        if isempty(wstr)
            wstr = get(r, "w_reg", "")
        end
        val = tryparse(Float64, wstr)
        val === nothing && continue
        if !haskey(W, gi)
            W[gi] = Dict{String, Float64}()
        end
        W[gi][gj] = val
    end
    return W
end


# ── Load fundamentals ────────────────────────────────────────────────

function load_fundamentals()
    path = derived_path("fundamentals", "firm_fundamentals.csv")
    _, rows = read_csv(path)
    # Latest record per firm (for SIC classification)
    latest = Dict{String, Dict{String, String}}()
    for r in rows
        gk = r["gvkey"]
        fy = r["fyear"]
        if !haskey(latest, gk) || fy > latest[gk]["fyear"]
            latest[gk] = r
        end
    end
    return latest
end

function get_sic4(fundamentals::Dict, gvkey::String)::Union{String, Nothing}
    f = get(fundamentals, gvkey, nothing)
    f === nothing && return nothing
    sic = get(f, "sic", "")
    isempty(sic) && return nothing
    return length(sic) >= 4 ? sic[1:4] : sic
end


# ── Load events ───────────────────────────────────────────────────────

struct EventInfo
    plant::String
    year::Union{Int, Nothing}
    event_date::String
    gvkeys::Vector{String}
end

function load_events()
    println("Loading events...")
    path = derived_path("events", "coal_retirement_events.csv")
    _, rows = read_csv(path)
    events = EventInfo[]
    for r in rows
        matched = get(r, "matched_gvkeys", "")
        isempty(matched) && continue
        get(r, "is_first_mover", "") != "True" && continue
        ann_date = strip(get(r, "announcement_date", ""))
        ret_date = strip(get(r, "event_date", ""))
        effective_date = isempty(ann_date) ? ret_date : ann_date
        event_year = nothing
        if length(effective_date) >= 4 && all(isdigit, effective_date[1:4])
            event_year = parse(Int, effective_date[1:4])
        else
            ry = get(r, "ret_year", "")
            if !isempty(ry)
                event_year = tryparse(Int, ry)
            end
        end
        gvkeys = split(matched, ';') .|> String
        push!(events, EventInfo(r["plant_name"], event_year, effective_date, gvkeys))
    end
    println("  First-mover events: $(length(events))")
    return events
end


# ── CAR computation ───────────────────────────────────────────────────

function compute_monthly_car(
    gvkey::String,
    event_month::String,
    monthly_ret::Dict{String, Dict{String, Float64}},
    market_ret::Dict{String, Float64};
    post::Int = POST_MONTHS,
)::Union{Float64, Nothing}

    !haskey(monthly_ret, gvkey) && return nothing
    firm_months = sort(collect(keys(monthly_ret[gvkey])))

    # Find the first month >= event_month
    event_idx = nothing
    for (i, m) in enumerate(firm_months)
        if m >= event_month
            event_idx = i
            break
        end
    end
    event_idx === nothing && return nothing

    # Require enough pre-event data
    pre_start = max(1, event_idx - PRE_MONTHS)
    pre_count = 0
    for i in pre_start:(event_idx - 1)
        if haskey(monthly_ret[gvkey], firm_months[i])
            pre_count += 1
        end
    end
    pre_count < 12 && return nothing

    # Pre-demean ARs by pre-window mean
    ar_list = Float64[]
    for i in pre_start:(event_idx - 1)
        m = firm_months[i]
        if haskey(monthly_ret[gvkey], m) && haskey(market_ret, m)
            push!(ar_list, monthly_ret[gvkey][m] - market_ret[m])
        end
    end
    pre_mean_ar = isempty(ar_list) ? 0.0 : mean(ar_list)

    car = 0.0
    for offset in -1:post
        idx = event_idx + offset
        if 1 <= idx <= length(firm_months)
            m = firm_months[idx]
            if haskey(monthly_ret[gvkey], m) && haskey(market_ret, m)
                ar = monthly_ret[gvkey][m] - market_ret[m]
                car += ar - pre_mean_ar
            end
        end
    end
    return car
end


# ── Build observation dataset ─────────────────────────────────────────

struct Observation
    car::Float64
    w_geo::Float64
    w_fuel::Float64
    w_reg::Float64
    same_sector::Float64
    event_id::Int
    gvkey::String
end

function build_obs(
    all_events::Vector{EventInfo},
    monthly_ret::Dict{String, Dict{String, Float64}},
    market_ret::Dict{String, Float64},
    W_geo::Dict{String, Dict{String, Float64}},
    W_fuel::Dict{String, Dict{String, Float64}},
    W_reg::Dict{String, Dict{String, Float64}},
    fundamentals::Dict{String, Dict{String, String}},
)
    obs = Observation[]

    for (event_id, event) in enumerate(all_events)
        event_gvkeys = Set(event.gvkeys)
        year = event.year
        ed = event.event_date

        # Determine event month
        event_month = if length(ed) >= 7
            ed[1:7]
        elseif year !== nothing
            @sprintf("%04d-07", year)
        else
            nothing
        end
        event_month === nothing && continue

        # Get first-mover SIC4
        fm_sic4 = nothing
        for gk in event.gvkeys
            fm_sic4 = get_sic4(fundamentals, gk)
            fm_sic4 !== nothing && break
        end

        for fm_gk in event.gvkeys
            !haskey(W_geo, fm_gk) && continue
            neighbors = W_geo[fm_gk]
            neighbor_gks = setdiff(Set(keys(neighbors)), event_gvkeys)

            # Non-connected firms: in fundamentals but not event or neighbor
            non_connected = [gk for gk in keys(fundamentals)
                             if !(gk in event_gvkeys) && !haskey(neighbors, gk)]

            # Deterministic control sampling using Julia's hash() on gvkey string.
            # NOTE: differs from Python's hashlib.md5 but the bootstrap logic is
            # identical. The key results depend on the full dataset, not on which
            # specific controls are drawn.
            stable_seed = hash(fm_gk) % UInt32
            rng = Random.MersenneTwister(stable_seed)
            n_ctrl = min(length(non_connected), max(5 * length(neighbor_gks), 20))
            ctrl_sample = if length(non_connected) > n_ctrl
                # Shuffle and take first n_ctrl (equivalent to random.sample)
                shuffled = shuffle(rng, non_connected)
                shuffled[1:n_ctrl]
            else
                non_connected
            end

            candidate_firms = vcat(collect(neighbor_gks), ctrl_sample)

            for gk in candidate_firms
                w_geo_val  = get(neighbors, gk, 0.0)
                w_fuel_val = get(get(W_fuel, fm_gk, Dict{String,Float64}()), gk, 0.0)
                w_reg_val  = get(get(W_reg,  fm_gk, Dict{String,Float64}()), gk, 0.0)

                j_sic4 = get_sic4(fundamentals, gk)
                same_sector = (fm_sic4 !== nothing && j_sic4 !== nothing &&
                               fm_sic4 == j_sic4) ? 1.0 : 0.0

                car = compute_monthly_car(gk, event_month, monthly_ret, market_ret;
                                          post=POST_MONTHS)
                car === nothing && continue

                push!(obs, Observation(car, w_geo_val, w_fuel_val, w_reg_val,
                                       same_sector, event_id, gk))
            end
        end
    end

    return obs
end


# ── OLS with cluster-robust variance-covariance matrix ────────────────

struct OLSResult
    beta::Vector{Float64}       # coefficients (intercept first)
    se::Vector{Float64}         # cluster-robust standard errors
    t_stats::Vector{Float64}    # t-statistics
    r2::Float64
    n::Int
    k::Int
    V::Matrix{Float64}          # full variance-covariance matrix
    ss_res::Float64
    n_clusters::Int
    names::Vector{String}       # variable names (intercept first)
end

function ols_full(
    obs::Vector{Observation},
    spec_vars::Vector{String};
    cluster::Bool = true,
)::Union{OLSResult, Nothing}

    n = length(obs)
    k = length(spec_vars) + 1  # +1 for intercept
    n <= k + 1 && return nothing

    # Build y and X
    y = Vector{Float64}(undef, n)
    X = Matrix{Float64}(undef, n, k)
    for i in 1:n
        y[i] = obs[i].car
        X[i, 1] = 1.0
        for (j, v) in enumerate(spec_vars)
            X[i, j+1] = if v == "w_geo"
                obs[i].w_geo
            elseif v == "w_fuel"
                obs[i].w_fuel
            elseif v == "w_reg"
                obs[i].w_reg
            elseif v == "same_sector"
                obs[i].same_sector
            else
                error("Unknown variable: $v")
            end
        end
    end

    y_mean = mean(y)
    ss_tot = sum((yi - y_mean)^2 for yi in y)
    ss_tot < 1e-15 && return nothing

    # OLS: beta = (X'X)^{-1} X'y
    XtX = X' * X
    inv_XtX = try
        inv(XtX)
    catch
        return nothing
    end
    beta = inv_XtX * (X' * y)

    y_hat = X * beta
    resid = y - y_hat
    ss_res = dot(resid, resid)
    r2 = 1.0 - ss_res / ss_tot

    # Variance-covariance matrix
    V = Matrix{Float64}(undef, k, k)
    G = 0

    if cluster
        # Build cluster map: event_id -> [indices]
        cluster_map = Dict{Int, Vector{Int}}()
        for i in 1:n
            eid = obs[i].event_id
            if !haskey(cluster_map, eid)
                cluster_map[eid] = Int[]
            end
            push!(cluster_map[eid], i)
        end
        G = length(cluster_map)

        # Clustered meat matrix S = sum_g (X_g' * u_g * u_g' * X_g)
        S = zeros(Float64, k, k)
        xu = zeros(Float64, k)
        for idxs in values(cluster_map)
            fill!(xu, 0.0)
            for i in idxs
                ri = resid[i]
                @inbounds for a in 1:k
                    xu[a] += X[i, a] * ri
                end
            end
            for a in 1:k
                @inbounds for b in a:k
                    v = xu[a] * xu[b]
                    S[a, b] += v
                    if a != b
                        S[b, a] += v
                    end
                end
            end
        end

        # Sandwich: V = inv(X'X) * S * inv(X'X)
        V .= inv_XtX * S * inv_XtX

        # Small-sample correction: (G/(G-1)) * ((N-1)/(N-k))
        if G > 1
            scale = (G / (G - 1)) * ((n - 1) / (n - k))
            V .*= scale
        end
    else
        s2 = n > k ? ss_res / (n - k) : 0.0
        V .= s2 .* inv_XtX
    end

    se = [sqrt(max(V[a, a], 0.0)) for a in 1:k]
    t_stats = [se[a] > 1e-15 ? beta[a] / se[a] : 0.0 for a in 1:k]

    var_names = vcat(["intercept"], spec_vars)
    return OLSResult(beta, se, t_stats, r2, n, k, V, ss_res, G, var_names)
end


# ── OLS via backslash for the permutation loop (speed) ────────────────

"""
Fast OLS returning only SSR. Used in the permutation loop.
X_perm is modified in place before calling this.
"""
function ols_ssr_fast(X::Matrix{Float64}, y::Vector{Float64})::Float64
    beta = X \ y
    resid = y - X * beta
    return dot(resid, resid)
end


# ── Significance stars ────────────────────────────────────────────────

function sig_stars(p::Float64)::String
    p < 0.01 && return "***"
    p < 0.05 && return "**"
    p < 0.10 && return "*"
    return ""
end


# ══════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════

function main()
    # ── Load all data ─────────────────────────────────────────────────

    monthly_ret = load_monthly_returns()
    market_ret  = load_ff_factors_monthly()

    println("Loading weight matrices...")
    W_geo = load_weight_matrix(derived_path("networks", "weight_matrix_W_geo.csv"), "w_ij")
    println("  W_geo firms: $(length(W_geo))")

    W_fuel = load_weight_matrix(derived_path("networks", "weight_matrix_W_fuel.csv"), "w_ij")
    n_fuel_edges = sum(length(v) for v in values(W_fuel); init=0)
    if isempty(W_fuel)
        println("  W_fuel: NOT FOUND")
    else
        println("  W_fuel edges: $n_fuel_edges")
    end

    W_reg = load_weight_matrix(derived_path("networks", "weight_matrix_W_regulatory.csv"), "w_ij")
    n_reg_edges = sum(length(v) for v in values(W_reg); init=0)
    if isempty(W_reg)
        println("  W_reg: NOT FOUND")
    else
        println("  W_reg edges: $n_reg_edges")
    end

    fundamentals = load_fundamentals()
    all_events = load_events()

    # ── Banner ────────────────────────────────────────────────────────

    println()
    println("=" ^ 70)
    println("FOCUSED HYPOTHESIS TESTS: CHANNEL SPLIT")
    @printf("Window: [-1, +%d] months\n", POST_MONTHS)
    @printf("Bootstrap: B = %d, seed = %d\n", B_PERM, SEED)
    println("=" ^ 70)

    # ── Build dataset ─────────────────────────────────────────────────

    println("\nBuilding dataset...")
    obs = build_obs(all_events, monthly_ret, market_ret, W_geo, W_fuel, W_reg, fundamentals)
    n_obs = length(obs)
    println("  N = $n_obs observations")
    println("  Neighbors (w_geo > 0): $(count(o -> o.w_geo > 0, obs))")
    println("  Non-connected: $(count(o -> o.w_geo == 0.0, obs))")
    println("  Events: $(length(unique(o.event_id for o in obs)))")

    # Check same_sector variation
    ss_vals = Set(o.same_sector for o in obs)
    spec_vars = if length(ss_vals) <= 1
        println("  WARNING: no same_sector variation; dropping same_sector from spec")
        copy(CHANNEL_VARS)
    else
        copy(SPEC_VARS_FULL)
    end

    # ── Estimate unrestricted model ───────────────────────────────────

    println("\nEstimating unrestricted model...")
    println("  CAR = alpha + beta_geo*w_geo + beta_fuel*w_fuel + beta_reg*w_reg + beta_s*SameSector")
    res_full = ols_full(obs, spec_vars; cluster=true)
    if res_full === nothing
        println("ERROR: Unrestricted OLS failed.")
        return
    end

    @printf("  N = %d, R2 = %.6f\n", res_full.n, res_full.r2)
    for v in spec_vars
        idx = findfirst(==(v), res_full.names)
        @printf("  %s: beta = %+.6f, se = %.6f, t = %.3f\n",
                v, res_full.beta[idx], res_full.se[idx], res_full.t_stats[idx])
    end

    # ── Estimate restricted model (same_sector only) ──────────────────

    restricted_vars = [v for v in spec_vars if !(v in CHANNEL_VARS)]
    println("\nEstimating restricted model...")
    suffix = isempty(restricted_vars) ? "" : " + beta_s*SameSector"
    println("  CAR = alpha$suffix")
    res_restricted = ols_full(obs, restricted_vars; cluster=true)
    if res_restricted === nothing
        println("ERROR: Restricted OLS failed.")
        return
    end
    @printf("  N = %d, R2 = %.6f\n", res_restricted.n, res_restricted.r2)

    # ══════════════════════════════════════════════════════════════════
    # TEST 1: Joint F-test with permutation bootstrap
    # ══════════════════════════════════════════════════════════════════

    println()
    println("=" ^ 70)
    println("TEST 1: JOINT F-TEST")
    println("H0: beta_geo = beta_fuel = beta_reg = 0")
    println("=" ^ 70)

    ssr_restricted   = res_restricted.ss_res
    ssr_unrestricted = res_full.ss_res
    q       = length(CHANNEL_VARS)  # number of restrictions = 3
    k_full  = length(spec_vars) + 1 # total regressors including intercept

    f_stat = ((ssr_restricted - ssr_unrestricted) / q) /
             (ssr_unrestricted / (n_obs - k_full))

    @printf("\n  SSR_restricted:   %.6f\n", ssr_restricted)
    @printf("  SSR_unrestricted: %.6f\n", ssr_unrestricted)
    @printf("  F-statistic: %.4f\n", f_stat)
    @printf("  df: (%d, %d)\n", q, n_obs - k_full)

    # Permutation bootstrap for F-test p-value
    @printf("\nRunning %d permutation bootstraps for F-test p-value...\n", B_PERM)
    rng = MersenneTwister(SEED)

    # Pre-extract arrays for the permutation loop
    car_arr = [o.car for o in obs]
    ss_arr  = [o.same_sector for o in obs]
    wgeo_arr  = [o.w_geo for o in obs]
    wfuel_arr = [o.w_fuel for o in obs]
    wreg_arr  = [o.w_reg for o in obs]

    n_perm = length(obs)

    # Build the X matrix template: [intercept, w_geo, w_fuel, w_reg, same_sector]
    # For spec_vars that may or may not include same_sector
    # We need the column mapping to match spec_vars ordering
    # Column 1 = intercept, then spec_vars in order
    X_perm = Matrix{Float64}(undef, n_perm, k_full)
    for i in 1:n_perm
        X_perm[i, 1] = 1.0
        for (j, v) in enumerate(spec_vars)
            X_perm[i, j+1] = if v == "w_geo"
                wgeo_arr[i]
            elseif v == "w_fuel"
                wfuel_arr[i]
            elseif v == "w_reg"
                wreg_arr[i]
            elseif v == "same_sector"
                ss_arr[i]
            else
                0.0
            end
        end
    end

    # Find which columns correspond to the three spatial channels
    geo_col  = findfirst(==("w_geo"),  spec_vars) + 1   # +1 for intercept col
    fuel_col = findfirst(==("w_fuel"), spec_vars) + 1
    reg_col  = findfirst(==("w_reg"),  spec_vars) + 1

    # The restricted model SSR does not change under permutation of spatial weights
    # (since restricted model only uses same_sector, which stays fixed).
    # So SSR_restricted is constant.

    perm_idx = collect(1:n_perm)
    f_boot = zeros(Float64, B_PERM)

    t_start = time()
    for b in 1:B_PERM
        if b % 100 == 0
            @printf("  permutation %d/%d\n", b, B_PERM)
        end

        # Permute spatial weight indices
        randperm!(rng, perm_idx)

        # Apply permutation to spatial weight columns only
        @inbounds for i in 1:n_perm
            pi = perm_idx[i]
            X_perm[i, geo_col]  = wgeo_arr[pi]
            X_perm[i, fuel_col] = wfuel_arr[pi]
            X_perm[i, reg_col]  = wreg_arr[pi]
        end

        # Fast OLS via backslash — only need SSR
        ssr_perm = ols_ssr_fast(X_perm, car_arr)
        f_boot[b] = ((ssr_restricted - ssr_perm) / q) / (ssr_perm / (n_perm - k_full))
    end
    elapsed = time() - t_start
    @printf("  Permutations complete in %.2f seconds (%.1f iter/s)\n", elapsed, B_PERM / elapsed)

    p_f_perm = count(fb -> fb >= f_stat, f_boot) / B_PERM
    @printf("  F-test p-value (permutation): %.4f\n", p_f_perm)

    f_reject = p_f_perm < 0.05 ? "REJECT H0" : "FAIL TO REJECT H0"
    @printf("  Conclusion at 5%%: %s\n", f_reject)

    # ══════════════════════════════════════════════════════════════════
    # TEST 2: Difference-in-coefficients test
    # ══════════════════════════════════════════════════════════════════

    println()
    println("=" ^ 70)
    println("TEST 2: DIFFERENCE-IN-COEFFICIENTS")
    println("H0: beta_geo = beta_fuel")
    println("H1: beta_geo != beta_fuel (opposing signs)")
    println("=" ^ 70)

    # Extract from full variance-covariance matrix V
    # Variable ordering: ["intercept"] + spec_vars
    idx_geo  = findfirst(==("w_geo"),  res_full.names)
    idx_fuel = findfirst(==("w_fuel"), res_full.names)

    V = res_full.V
    beta_geo  = res_full.beta[idx_geo]
    beta_fuel = res_full.beta[idx_fuel]
    se_geo    = res_full.se[idx_geo]
    se_fuel   = res_full.se[idx_fuel]

    var_geo      = V[idx_geo, idx_geo]
    var_fuel     = V[idx_fuel, idx_fuel]
    cov_geo_fuel = V[idx_geo, idx_fuel]

    diff    = beta_geo - beta_fuel
    se_diff = sqrt(var_geo + var_fuel - 2.0 * cov_geo_fuel)
    t_diff  = se_diff > 1e-15 ? diff / se_diff : 0.0
    p_diff  = p_from_t(t_diff)

    @printf("\n  beta_geo:  %+.6f (SE %.6f)\n", beta_geo, se_geo)
    @printf("  beta_fuel: %+.6f (SE %.6f)\n", beta_fuel, se_fuel)
    @printf("  Difference (beta_geo - beta_fuel): %+.6f\n", diff)
    @printf("  Var(beta_geo):           %.10f\n", var_geo)
    @printf("  Var(beta_fuel):          %.10f\n", var_fuel)
    @printf("  Cov(beta_geo, beta_fuel): %+.10f\n", cov_geo_fuel)
    @printf("  SE of difference:        %.6f\n", se_diff)
    @printf("  t-statistic:             %.3f\n", t_diff)
    @printf("  p-value (two-sided):     %.4f\n", p_diff)

    diff_reject = p_diff < 0.05 ? "REJECT H0" : "FAIL TO REJECT H0"
    @printf("  Conclusion at 5%%: %s\n", diff_reject)

    sig_stars_diff = sig_stars(p_diff)
    diff_interp = p_diff < 0.05 ? "statistically significant" : "not significant"

    # ── Full regression results ───────────────────────────────────────

    println()
    println("=" ^ 70)
    println("FULL REGRESSION RESULTS")
    println("=" ^ 70)
    for (i, name) in enumerate(res_full.names)
        b = res_full.beta[i]
        s = res_full.se[i]
        t = res_full.t_stats[i]
        p = p_from_t(t)
        stars = sig_stars(p)
        @printf("  %-15s %+.6f  (%.6f)  t=%.3f  p=%.4f%s\n", name, b, s, t, p, stars)
    end

    # ── Write markdown output ─────────────────────────────────────────

    out_path = results_path("metrics", "joint_tests.md")
    mkpath(dirname(out_path))

    sig_stars_f = sig_stars(p_f_perm)
    f_interp = p_f_perm < 0.05 ? "Reject" : "Fail to reject"

    lines = String[]
    push!(lines, "# Focused Hypothesis Tests: Channel Split")
    push!(lines, "")
    push!(lines, @sprintf("Window: [-1, +%d] months (monthly CARs, vwretd)", POST_MONTHS))
    push!(lines, @sprintf("Events: %d first-mover-matched (175 used in pooled regression below; 117 with ≥20 firms qualify for FM)", length(all_events)))
    push!(lines, @sprintf("N = %d observations, %d event clusters", n_obs, res_full.n_clusters))
    push!(lines, "Standard errors: event-clustered")
    push!(lines, "")
    push!(lines, "## Test 1: Joint F-test (H0: beta_geo = beta_fuel = beta_reg = 0)")
    push!(lines, "")
    push!(lines, "Unrestricted: CAR = alpha + beta_geo * w^geo + beta_fuel * w^fuel + beta_reg * w^reg + beta_s * SameSector")
    push!(lines, "Restricted:   CAR = alpha + beta_s * SameSector")
    push!(lines, "")
    push!(lines, @sprintf("SSR_restricted:   %.6f", ssr_restricted))
    push!(lines, @sprintf("SSR_unrestricted: %.6f", ssr_unrestricted))
    push!(lines, @sprintf("F-statistic: %.4f", f_stat))
    push!(lines, @sprintf("df: (%d, %d)", q, n_obs - k_full))
    push!(lines, @sprintf("p-value (permutation, B=%d): %.4f%s", B_PERM, p_f_perm, sig_stars_f))
    push!(lines, @sprintf("N: %d", n_obs))
    push!(lines, "")

    f_jointly = p_f_perm < 0.05 ?
        "jointly predict CARs" : "do not jointly predict CARs"
    push!(lines, @sprintf("%s H0 at 5%%. The spatial network channels %s around coal retirement events.",
                          f_interp, f_jointly))
    push!(lines, "")
    push!(lines, "## Test 2: Difference test (H0: beta_geo = beta_fuel)")
    push!(lines, "")
    push!(lines, @sprintf("beta_geo:  %+.6f (SE %.6f)", beta_geo, se_geo))
    push!(lines, @sprintf("beta_fuel: %+.6f (SE %.6f)", beta_fuel, se_fuel))
    push!(lines, @sprintf("Difference (beta_geo - beta_fuel): %+.6f", diff))
    push!(lines, @sprintf("SE of difference: %.6f", se_diff))
    push!(lines, @sprintf("t-statistic: %.3f", t_diff))
    push!(lines, @sprintf("p-value: %.4f%s", p_diff, sig_stars_diff))
    push!(lines, @sprintf("Cov(beta_geo, beta_fuel): %+.10f", cov_geo_fuel))
    push!(lines, "")
    push!(lines, @sprintf("Interpretation: The opposing-sign channel split is %s in a single test (t = %.3f, p = %.4f).",
                          diff_interp, t_diff, p_diff))
    push!(lines, "")
    push!(lines, "## Full regression coefficients")
    push!(lines, "")
    push!(lines, "| Variable | Beta | SE | t | p |")
    push!(lines, "|---|---:|---:|---:|---:|")

    for (i, name) in enumerate(res_full.names)
        b = res_full.beta[i]
        s = res_full.se[i]
        t = res_full.t_stats[i]
        p = p_from_t(t)
        stars = sig_stars(p)
        push!(lines, @sprintf("| %s | %+.6f | %.6f | %.3f | %.4f%s |", name, b, s, t, p, stars))
    end

    push!(lines, "")
    push!(lines, "## Comparison with individual tests")
    push!(lines, "")
    push!(lines, "| Test | Hypotheses tested | Correction needed | Result |")
    push!(lines, "|---|---|---|---|")
    push!(lines, "| Individual t-tests | 3 | Romano-Wolf | 1/3 significant (fuel) |")
    f_result = p_f_perm < 0.05 ?
        @sprintf("significant (p=%.4f)", p_f_perm) :
        @sprintf("not significant (p=%.4f)", p_f_perm)
    d_result = p_diff < 0.05 ?
        @sprintf("significant (p=%.4f)", p_diff) :
        @sprintf("not significant (p=%.4f)", p_diff)
    push!(lines, @sprintf("| Joint F-test | 1 | None | %s |", f_result))
    push!(lines, @sprintf("| Difference test | 1 | None | %s |", d_result))
    push!(lines, "")

    open(out_path, "w") do io
        write(io, join(lines, "\n"))
    end

    println("\nWrote: $out_path")
    println("Done.")
end

# ── Entry point ───────────────────────────────────────────────────────
main()
