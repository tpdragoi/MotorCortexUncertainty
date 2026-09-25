# Library Import
using Random
using Distributions
using Plots
using StatsBase
using CSV
using DataFrames
using LinearAlgebra
using StateSpaceDynamics
using Statistics
const SSD = StateSpaceDynamics

# Column order of the Keypoint_Feats_*.csv files, matching traj_features in
# MATLAB/Preprocess_Data.m.
const KP_FEATURE_NAMES = [
    "tongue_xdisp_view1", "left_tongue_xdisp_view1", "right_tongue_xdisp_view1",
    "jaw_xdisp_view1", "trident_xdisp_view1", "nose_xdisp_view1", "tongue_ydisp_view1",
    "left_tongue_ydisp_view1", "right_tongue_ydisp_view1", "jaw_ydisp_view1", "trident_ydisp_view1",
    "nose_ydisp_view1", "top_tongue_xdisp_view2", "topleft_tongue_xdisp_view2",
    "bottom_tongue_xdisp_view2", "bottomleft_tongue_xdisp_view2", "jaw_xdisp_view2",
    "top_nostril_xdisp_view2", "bottom_nostril_xdisp_view2", "top_tongue_ydisp_view2",
    "topleft_tongue_ydisp_view2", "bottom_tongue_ydisp_view2", "bottomleft_tongue_ydisp_view2",
    "jaw_ydisp_view2", "top_nostril_ydisp_view2", "bottom_nostril_ydisp_view2",
]

# Time axis (seconds, relative to true GC) for every trial's 600-point Uncut
# chunk from Preprocess_Data.m: tmin=-2.0, tmax=5.0, dt=1/100, then the
# `101:end` trim
const TIME_UNCUT = collect(range(-0.9950, 4.9950; length=600))

"""
Define custom fit functions for ridge regression using StateSpaceDynamics.jl
The default models set the ridge parameter to 0.0
"""
function fit_custom!(
    model::HiddenMarkovModel,
    Y::Vector{<:Matrix{<:Real}},
    X::Union{Vector{<:Matrix{<:Real}},Nothing}=nothing;
    max_iters::Int=100,
    tol::Float64=1e-6,
)
    lls = [-Inf]
    data = X === nothing ? (Y,) : (X, Y)

    # Initialize log_likelihood
    log_likelihood = -Inf

    # Transform each matrix in each tuple to the correct orientation
    transposed_matrices = map(data_tuple -> Matrix.(transpose.(data_tuple)), data)
    zipped_matrices = collect(zip(transposed_matrices...))
    total_obs = sum(size(trial_mat[1], 1) for trial_mat in zipped_matrices)

    # initialize a vector of ForwardBackward storage and an aggregate storage
    FB_storage_vec = [SSD.initialize_forward_backward(model, size(trial_tuple[1],1)) for trial_tuple in zipped_matrices]
    Aggregate_FB_storage = SSD.initialize_forward_backward(model, total_obs)

    for iter in 1:max_iters
        println("A: ", model.A)
        println("Iter: ", iter)
        # broadcast estep!() to all storage structs
        output = SSD.estep!.(Ref(model), zipped_matrices, FB_storage_vec)

        # collect storage stucts into one struct for m step
        SSD.aggregate_forward_backward!(Aggregate_FB_storage, FB_storage_vec)

        # Calculate log_likelihood
        log_likelihood_current = sum([SSD.logsumexp(FB_vec.α[:, end]) / size(FB_vec.α, 2) for FB_vec in FB_storage_vec])
        push!(lls, log_likelihood_current)

        # Check for convergence
        if abs(log_likelihood_current - log_likelihood) < tol
            break
        else
            log_likelihood = log_likelihood_current
        end

        # Get data trial tuples stacked for mstep!()
        stacked_data = SSD.stack_tuples(zipped_matrices)

        # M_step
        SSD.mstep!(model, FB_storage_vec, Aggregate_FB_storage, stacked_data)
    end

    return lls
end

function SSD.update_emissions!(model::SSD.AbstractHMM, FB_storage::SSD.ForwardBackward, data)
    # Use closed-form weighted ridge regression instead of SSD's default optimizer
    w = exp.(permutedims(FB_storage.γ))

    for k in 1:(model.K)
        β, Σ = weighted_ridge_regression(data..., model.B[1].λ, w=w[:,k])
        model.B[k].β = β
        model.B[k].Σ = Σ
    end
end

function weighted_ridge_regression(
    X::Matrix{Float64},
    Y::Matrix{Float64},
    λ::Float64;
    w::Vector{Float64}=ones(size(X, 1))
)
    @assert size(X, 1) == size(Y, 1) == length(w)

    N, D = size(X)
    X_bias = hcat(ones(N), X)
    W = Diagonal(w)

    A = X_bias' * W * X_bias
    b = X_bias' * W * Y

    # Ridge regularization (don't penalize intercept)
    reg = zeros(D + 1, D + 1)
    reg[2:end, 2:end] .= λ

    β = (A + reg) \ b

    residuals = Y - X_bias * β
    Σ = (residuals' * W * residuals) / sum(w)
    Σ = 0.5 * (Σ + Σ')  # Force symmetry

    return β, Σ
end

function r2_score(y_true, y_pred)
    ss_res = sum((y_true .- y_pred).^2)
    ss_tot = sum((y_true .- mean(y_true, dims=1)).^2)
    return 1 - (ss_res / ss_tot)
end

function label_data(model, Y::Vector{<:Matrix{<:Real}}, X::Union{Vector{<:Matrix{<:Real}}, Nothing}=nothing)
    data = X === nothing ? (Y,) : (X, Y)

    # Transpose each matrix to shape (features × timepoints)
    transposed_matrices = map(data_tuple -> Matrix.(transpose.(data_tuple)), data)
    zipped_matrices = collect(zip(transposed_matrices...))  # Vector of tuples: (Xᵀ, Yᵀ)

    total_obs = sum(size(trial_tuple[1], 1) for trial_tuple in zipped_matrices)

    # Allocate ForwardBackward storage
    FB_storage_vec = [SSD.initialize_forward_backward(model, size(trial_tuple[1], 1)) for trial_tuple in zipped_matrices]
    Aggregate_FB_storage = SSD.initialize_forward_backward(model, total_obs)

    # Run E-step
    SSD.estep!.(Ref(model), zipped_matrices, FB_storage_vec)

    return FB_storage_vec
end

"""
Data loading functions for full trials for learning and chopped trials for inference.
"""
function load_data(path, condition, probenum, chunk_size)
    # Helper function to chunk matrix into 600-timepoint segments
    function chunk_matrix(mat, chunk_size)
        num_chunks = size(mat, 1) ÷ chunk_size
        [mat[(i-1)*chunk_size + 1 : i*chunk_size, :] for i in 1:num_chunks]
    end

    PCA_P1_path = path * "PCA_Probe" * string(probenum) * "_" * condition * "_Uncut.csv"
    KP_path = path * "Keypoint_Feats_" * condition * "_Uncut.csv"

    # In the previous step if the session only used probe2, it was saved to probe 1 by default.
    # This is detected if the probenum is 2 but the PCA for probe 2 doesn't exist.
    if probenum == 2 && !isfile(PCA_P1_path)
        PCA_P1_path = path * "PCA_Probe1_" * condition * "_Uncut.csv"
    end

    # Load the data into matrices
    PCA_P1_mat = Matrix(CSV.read(PCA_P1_path, DataFrame; header=false))
    KP_mat = Matrix(CSV.read(KP_path, DataFrame; header=false))

    # Chunk all matrices
    PCA_P1_chunks = chunk_matrix(PCA_P1_mat, chunk_size)
    KP_chunks = chunk_matrix(KP_mat, chunk_size)

    return PCA_P1_chunks, KP_chunks
end

function load_data_cut(path, condition, probe_num)
    # Construct file paths
    PCA_P1_path = path * "PCA_Probe" * string(probe_num) * "_" * condition * "_Cut.csv"
    KP_path = path * "Keypoint_Feats_" * condition * "_Cut.csv"

    # In the previous step if the session only used probe2, it was saved to probe 1 by default.
    # This is detected if the probenum is 2 but the PCA for probe 2 doesn't exist.
    if probe_num == 2 && !isfile(PCA_P1_path)
        PCA_P1_path = path * "PCA_Probe1_" * condition * "_Cut.csv"
    end
    FCs_path = path * "FCs_" * condition * ".csv"
    LRCs_path = path * "LRCs_" * condition * ".csv"
    Tongue_path = path * "Tongue_" * condition * ".csv"

    # Load the data into matrices
    PCA_P1_mat = Matrix(CSV.read(PCA_P1_path, DataFrame; header=false))
    KP_mat = Matrix(CSV.read(KP_path, DataFrame; header=false))
    FCs_mat = Matrix(CSV.read(FCs_path, DataFrame; header=false))
    LRCs = vec(Matrix(CSV.read(LRCs_path, DataFrame; header=false)))
    Tongue_mat = Matrix(CSV.read(Tongue_path, DataFrame; header=false))

    # Chunking function
    function chunk_matrix(mat, lengths)
        chunks = Vector{Matrix{Float64}}()
        start_idx = 1
        for len in lengths
            stop_idx = start_idx + len - 1
            push!(chunks, mat[start_idx:stop_idx, :])
            start_idx = stop_idx + 1
        end
        return chunks
    end

    # Chunk all data at SR Hz
    PCA_P1_chunks = chunk_matrix(PCA_P1_mat, LRCs)
    KP_chunks = chunk_matrix(KP_mat, LRCs)

    return PCA_P1_chunks, KP_chunks, FCs_mat, LRCs, Tongue_mat
end


"""
Functions for prepping and making design matrices for the fitting procedure.
"""
function kernelize_window_features(X_train::Vector{Matrix{T}}, lags::Int=4, leads::Int=0) where T
    processed_features = Vector{Matrix{T}}(undef, length(X_train))

    for (i, X) in enumerate(X_train)
        num_timepoints, num_features = size(X)

        # Define valid time range where full lag and lead context exists
        start_idx = lags + 1
        end_idx = num_timepoints - leads

        # Initialize new feature matrix
        num_new_timepoints = end_idx - start_idx + 1
        windowed_features = Matrix{T}(undef, num_new_timepoints, num_features * (lags + leads + 1))

        for (j, t) in enumerate(start_idx:end_idx)
            # Collect lagged, current, and leading feature values
            feature_vector = vcat(X[t, :],
                                    (X[t - lag, :] for lag in lags:-1:1)...,   # past (lags -> 1)
                                                                  # current
                                  (X[t + lead, :] for lead in 1:leads)...)    # future (1 -> leads)
            windowed_features[j, :] = feature_vector
        end

        processed_features[i] = windowed_features
    end

    return processed_features
end

function detect_probe(session::String)
    endswith(session, "_P3") && return 3
    endswith(session, "_P1") && return 1
    endswith(session, "_P2") && return 2
    return nothing
end

function prep_session_data(session_path::String, prb::Int, sr::Int)
    PCA_R1, KP_R1 = load_data(session_path, "R1", prb, sr)
    PCA_R4, KP_R4 = load_data(session_path, "R4", prb, sr)

    PCA_R1_Cut, KP_R1_Cut, FCs_R1, _, Tongue_R1 = load_data_cut(session_path, "R1", prb)
    PCA_R4_Cut, KP_R4_Cut, FCs_R4, _, Tongue_R4 = load_data_cut(session_path, "R4", prb)

    @assert length(TIME_UNCUT) == sr "TIME_UNCUT has $(length(TIME_UNCUT)) points but CHUNK=$sr -- update TIME_UNCUT if the chunk size changes"
    time_vec = TIME_UNCUT

    # Replace lingering NaN values with zeros
    for KP in (KP_R1, KP_R4, KP_R1_Cut, KP_R4_Cut)
        fill_nan_zero!(KP)
    end

    return (; PCA_R1, KP_R1, PCA_R4, KP_R4,
             PCA_R1_Cut, KP_R1_Cut, FCs_R1,
             PCA_R4_Cut, KP_R4_Cut, FCs_R4,
             Tongue_R1, Tongue_R4, time_vec)
end

# Replace NaN with 0.0 in every trial matrix in a vector of KP matrices, in place.
function fill_nan_zero!(mats::Vector{<:Matrix})
    for i in eachindex(mats)
        mats[i] .= replace(mats[i], NaN => 0.0)
    end
    return mats
end

# Drop the given column indices from every trial matrix in a vector of KP matrices, in place.
function drop_columns!(mats::Vector{<:Matrix}, cols_to_drop)
    keep = setdiff(1:size(mats[1], 2), cols_to_drop)
    for i in eachindex(mats)
        mats[i] = mats[i][:, keep]
    end
    return mats
end

# Prefit the engaged encoder to initialize the GLM-HMM fitting procedure.
function prefit_engaged_encoder(KP_R1, KP_R4, PCA_R1, PCA_R4, FCs_R1, FCs_R4;
                                 lags, start_time, pre, post)
    X_R1 = [X[start_time-lags:end, :] for X in KP_R1]
    X_R4 = [X[start_time-lags:end, :] for X in KP_R4]

    Y_R1 = [Y[start_time-lags:end, 1:10] for Y in PCA_R1]
    Y_R4 = [Y[start_time-lags:end, 1:10] for Y in PCA_R4]

    X_R1_kernel = kernelize_window_features(X_R1, lags)
    X_R4_kernel = kernelize_window_features(X_R4, lags)

    Y_R1_trimmed = kernelize_window_features(Y_R1, lags)
    Y_R4_trimmed = kernelize_window_features(Y_R4, lags)

    FCs_R1 = FCs_R1 .- start_time
    FCs_R4 = FCs_R4 .- start_time

    X_R1 = [X_R1_kernel[i][(FCs_R1[i]-pre):(FCs_R1[i]), :] for i in eachindex(X_R1_kernel)]
    X_R4 = [X_R4_kernel[i][(FCs_R4[i]-pre):(FCs_R4[i]+post), :] for i in eachindex(X_R4_kernel)]

    Y_R1 = [Y_R1_trimmed[i][(FCs_R1[i]-pre):(FCs_R1[i]), :] for i in eachindex(Y_R1_trimmed)]
    Y_R4 = [Y_R4_trimmed[i][(FCs_R4[i]-pre):(FCs_R4[i]+post), :] for i in eachindex(Y_R4_trimmed)]

    X_eng = cat(X_R1, X_R4, dims=1)
    Y_eng = cat(Y_R1, Y_R4, dims=1)

    X_eng = vcat(X_eng...)
    Y_eng = vcat(Y_eng...)

    λ_ridge = 0.0

    β_eng, Σ_eng = weighted_ridge_regression(X_eng, Y_eng, λ_ridge)
    return β_eng, Σ_eng
end

"""Build the switching-model design matrices from the FC-cut trials and fit the
2-state Gaussian HMM-GLM, seeding state 1 with the prefit engaged encoder"""
function fit_switching_model(KP_R1_Cut, KP_R4_Cut, PCA_R1_Cut, PCA_R4_Cut, β_eng, Σ_eng;
                              lags, start_time, max_iters)
    X_R1 = [X[start_time-lags:end, :] for X in KP_R1_Cut]
    X_R4 = [X[start_time-lags:end, :] for X in KP_R4_Cut]
    X = cat(X_R1, X_R4, dims=1)

    Y = cat(PCA_R1_Cut, PCA_R4_Cut, dims=1)
    Y = [y[start_time-lags:end, 1:10] for y in Y]

    X_kern = kernelize_window_features(X, lags)
    Y_trim = kernelize_window_features(Y, lags)

    X_ready = permutedims.(X_kern)
    Y_ready = permutedims.(Y_trim)

    model = SwitchingGaussianRegression(;
        K=2, input_dim=size(X_ready[1])[1], output_dim=size(Y_ready[1])[1], include_intercept=true)

    model.B[1].β = β_eng
    model.B[1].Σ = Σ_eng
    model.B[1].λ = 0.0
    model.B[2].λ = 0.0

    model.A = [0.99 0.01; 0.01 0.99]
    model.πₖ = [0.01; 0.99]

    lls = fit_custom!(model, Y_ready, X_ready, max_iters=max_iters)

    return model, lls
end

"""Run forward-backward & Viterbi over each uncut trial (through `three_sec`) to
get per-timepoint engaged-state probabilities and hard state labels."""
function compute_state_inference(model, KP_R1, KP_R4, PCA_R1, PCA_R4, Tongue_R1, Tongue_R4, time_vec;
                                  lags, start_time, three_sec, two_sec)

    time_reg = time_vec[start_time:three_sec]

    X_R1 = [X[start_time-lags:three_sec, :] for X in KP_R1]
    X_R4 = [X[start_time-lags:three_sec, :] for X in KP_R4]

    Y_R1 = [Y[start_time-lags:three_sec, 1:10] for Y in PCA_R1]
    Y_R4 = [Y[start_time-lags:three_sec, 1:10] for Y in PCA_R4]

    X_R1_kernel = kernelize_window_features(X_R1, lags)
    X_R4_kernel = kernelize_window_features(X_R4, lags)

    Y_R1_trimmed = kernelize_window_features(Y_R1, lags)
    Y_R4_trimmed = kernelize_window_features(Y_R4, lags)

    YY    = permutedims.(Y_R1_trimmed)
    XX    = permutedims.(X_R1_kernel)
    YY_R4 = permutedims.(Y_R4_trimmed)
    XX_R4 = permutedims.(X_R4_kernel)

    FB_R1 = label_data(model, YY, XX)
    FB_R4 = label_data(model, YY_R4, XX_R4)

    V1 = SSD.viterbi(model, YY, XX)
    V4 = SSD.viterbi(model, YY_R4, XX_R4)

    γ_vectors_R1 = [FB_R1[K].γ[1, :] for K in eachindex(FB_R1)]
    γ_mean_R1 = mean(exp.(hcat(γ_vectors_R1...)), dims=2)

    γ_vectors_R4 = [FB_R4[K].γ[1, :] for K in eachindex(FB_R4)]
    γ_mean_R4 = mean(exp.(hcat(γ_vectors_R4...)), dims=2)

    Tongue_R1_win = Tongue_R1[start_time:three_sec, :]
    Tongue_R4_win = Tongue_R4[start_time:three_sec, :]

    R4_States = permutedims(hcat(γ_vectors_R4...))
    R1_States = permutedims(hcat(γ_vectors_R1...))

    R4_Vit = permutedims(hcat(V4...))
    R1_Vit = permutedims(hcat(V1...))

    # Plus 10 here to include the 10 bins pre GC. This will make this 211 points total
    R4_Tongue_df = DataFrame(permutedims(Tongue_R4_win[1:two_sec+11, :]), :auto)
    R1_Tongue_df = DataFrame(permutedims(Tongue_R1_win[1:two_sec+11, :]), :auto)

    return (; X_R1_kernel, X_R4_kernel, Y_R1_trimmed, Y_R4_trimmed,
             γ_mean_R1, γ_mean_R4,
             R1_States, R4_States, R1_Vit, R4_Vit,
             R1_Tongue_df, R4_Tongue_df, time_reg)
end

"""Per-PC R² and NMSE broken out by engaged vs. disengaged state.
Prints a summary table for `session`."""
function compute_encoding_accuracy(model, X_R1_kernel, X_R4_kernel, Y_R1_trimmed, Y_R4_trimmed,
                                    R1_States, R4_States, session; n_pcs)
    _, O = size(Y_R4_trimmed[1])

    ss_err_eng = zeros(O);  n_eng  = zeros(Int, O)
    ss_err_dis = zeros(O);  n_dis  = zeros(Int, O)
    sum_y = zeros(O);       sum_y2 = zeros(O);  n_tot = 0

    function accumulate!(X_kernel, Y_trimmed, States, r2_scores)
        for trial in eachindex(X_kernel)
            X_trial = X_kernel[trial]
            Y_trial = Y_trimmed[trial]
            T = size(X_trial, 1)
            y_pred = zeros(T, O)
            X_bias = hcat(ones(T), X_trial)
            for i in 1:T
                state = exp(States[trial, i])
                if state == 1.0
                    y_pred[i, :] = reshape(X_bias[i, :], 1, :) * model.B[1].β
                    for d in 1:O
                        ss_err_eng[d] += (Y_trial[i,d] - y_pred[i,d])^2
                        n_eng[d] += 1
                    end
                else
                    y_pred[i, :] = reshape(X_bias[i, :], 1, :) * model.B[2].β
                    for d in 1:O
                        ss_err_dis[d] += (Y_trial[i,d] - y_pred[i,d])^2
                        n_dis[d] += 1
                    end
                end
                for d in 1:O
                    sum_y[d]  += Y_trial[i,d]
                    sum_y2[d] += Y_trial[i,d]^2
                end
                n_tot += 1
            end
            for pc in 1:O
                r2_scores[trial, pc] = r2_score(Y_trial[:, pc], y_pred[:, pc])
            end
        end
    end

    r2_scores_R4 = zeros(length(X_R4_kernel), O)
    accumulate!(X_R4_kernel, Y_R4_trimmed, R4_States, r2_scores_R4)

    r2_scores_R1 = zeros(length(X_R1_kernel), O)
    accumulate!(X_R1_kernel, Y_R1_trimmed, R1_States, r2_scores_R1)

    # Overall R² per output dim (mean/std across trials)
    r2_all = vcat(r2_scores_R1, r2_scores_R4)
    mean_r2_per_pc = vec(mean(r2_all, dims=1))
    std_r2_per_pc  = vec(std(r2_all, dims=1))

    # NMSE = MSE_state / Var(y_overall) — shared denominator makes states comparable
    var_y   = sum_y2 ./ n_tot .- (sum_y ./ n_tot).^2
    nmse_eng = (ss_err_eng ./ max.(n_eng, 1)) ./ var_y
    nmse_dis = (ss_err_dis ./ max.(n_dis, 1)) ./ var_y

    # Average over lag dims to get one value per original PC
    n_lags_pc = O ÷ n_pcs
    mean_r2_pc  = [mean([mean_r2_per_pc[pc + n_pcs * lag] for lag in 0:(n_lags_pc-1)]) for pc in 1:n_pcs]
    std_r2_pc   = [mean([std_r2_per_pc[pc  + n_pcs * lag] for lag in 0:(n_lags_pc-1)]) for pc in 1:n_pcs]
    nmse_eng_pc = [mean([nmse_eng[pc + n_pcs * lag] for lag in 0:(n_lags_pc-1)]) for pc in 1:n_pcs]
    nmse_dis_pc = [mean([nmse_dis[pc + n_pcs * lag] for lag in 0:(n_lags_pc-1)]) for pc in 1:n_pcs]

    println("\n--- Encoding Accuracy (R²) per PC for session: $session ---")
    println("       " * rpad("Mean R²", 12) * "Std R²")
    println("       " * repeat("-", 24))
    for pc in 1:n_pcs
        println("PC $(lpad(pc, 2)): $(rpad(round(mean_r2_pc[pc], digits=3), 12)) $(round(std_r2_pc[pc], digits=3))")
    end

    println("\n--- NMSE by State ---")
    println("       " * rpad("Engaged", 12) * "Disengaged")
    println("       " * repeat("-", 24))
    for pc in 1:n_pcs
        println("PC $(lpad(pc, 2)): $(rpad(round(nmse_eng_pc[pc], digits=3), 12)) $(round(nmse_dis_pc[pc], digits=3))")
    end
    println("---------------------------------------------------\n")

    return (; mean_r2_per_pc, mean_r2_pc, std_r2_pc, nmse_eng_pc, nmse_dis_pc)
end

"""Plot the R²/NMSE summary and write it + the per-trial state/tongue data
to `output_dir/session_save/`."""
function save_session_results(output_dir, session_save, session, accuracy, inference; n_pcs)
    (; mean_r2_per_pc, mean_r2_pc, std_r2_pc, nmse_eng_pc, nmse_dis_pc) = accuracy
    (; R1_States, R4_States, R1_Vit, R4_Vit, R1_Tongue_df, R4_Tongue_df, time_reg) = inference

    session_dir = joinpath(output_dir, session_save)
    if !isdir(session_dir)
        mkpath(session_dir)
    end

    p_r2 = bar(1:n_pcs, mean_r2_pc;
        yerr        = std_r2_pc,
        xlabel      = "PC",
        ylabel      = "R²",
        title       = "Encoding R² per PC\n$session",
        label       = "Mean ± Std",
        ylims       = (0, 1),
        xticks      = 1:n_pcs,
        color       = :steelblue,
        alpha       = 0.8)

    nmse_ymax = max(maximum(nmse_eng_pc), maximum(nmse_dis_pc)) * 1.15
    p_nmse = plot(1:n_pcs, nmse_eng_pc;
        marker      = :circle,
        label       = "Engaged",
        xlabel      = "PC",
        ylabel      = "NMSE",
        title       = "NMSE by State\n$session",
        ylims       = (0, nmse_ymax),
        xticks      = 1:n_pcs)
    plot!(1:n_pcs, nmse_dis_pc; marker=:square, label="Disengaged")
    hline!([1.0]; linestyle=:dash, color=:black, label="Chance (NMSE=1)")

    p_combined = plot(p_r2, p_nmse; layout=(1, 2), size=(1000, 420))
    display(p_combined)
    savefig(p_combined, joinpath(session_dir, "encoding_accuracy.png"))

    println("**SAVE PATH**", joinpath(session_dir, "R14_PC_R2_Reg.csv"))

    R4_States_Vit_df = DataFrame(R4_Vit, :auto)
    R1_States_Vit_df = DataFrame(R1_Vit, :auto)
    R4_States_df = DataFrame(R4_States, :auto)
    R1_States_df = DataFrame(R1_States, :auto)
    mean_r2_df = DataFrame(mean_r2_per_pc', :auto)  # make it a 1×n_pcs DataFrame

    CSV.write(joinpath(session_dir, "R14_PC_R2_Reg.csv"), mean_r2_df; header=false)
    CSV.write(joinpath(session_dir, "R14_States_Reg.csv"), R4_States_df; header=false)
    CSV.write(joinpath(session_dir, "R1_States_Reg.csv"), R1_States_df; header=false)
    CSV.write(joinpath(session_dir, "R14_States_Vit_Reg.csv"), R4_States_Vit_df; header=false)
    CSV.write(joinpath(session_dir, "R1_States_Vit_Reg.csv"), R1_States_Vit_df; header=false)
    CSV.write(joinpath(session_dir, "R14_Tongue_Reg.csv"), R4_Tongue_df; header=false)
    CSV.write(joinpath(session_dir, "R1_Tongue_Reg.csv"), R1_Tongue_df; header=false)

    # Seconds, relative to true GC, one value per column of R1_States_Reg.csv /
    # R14_States_Reg.csv (shared by both -- see compute_state_inference).
    CSV.write(joinpath(session_dir, "Time_Reg.csv"), DataFrame(time = time_reg); header=false)

    println("SESSION DATA SAVED")
end
