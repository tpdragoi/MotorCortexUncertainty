# Library and Utilities Import
using Glob
include(".\\zutils.jl")

# Configuration (Setup paths according to local machine)
const BASE_PATH  = ".\\Processed_Sessions"
const OUTPUT_DIR = ".\\Results"

# Windowing / training parameters currently configured for 10ms bins.
const CHUNK      = 600  # chunk size
const LAGS       = 4
const START_TIME = 90
const PRE        = 3  # 3 bins pre GC for initialization
const POST       = 10  # 10 bins post GC for initialization
const THREE_SEC  = 300
const TWO_SEC    = 200
const N_PCS      = 10
const MAX_ITERS  = 100
const NOSTRIL_COLS = findall(n -> occursin("nostril", n), KP_FEATURE_NAMES)

# Loop over the sessions and fit the GLM-HMM model to each one independently
session_folders = filter(isdir, glob("*", BASE_PATH))

for session_path in session_folders
    session = splitpath(session_path)[end]
    session_save = replace(session, "-" => "_")
    println("session save: ", session_save)

    # Skip folders without probe labels
    prb = detect_probe(session)
    if prb === nothing
        println("Skipping folder without probe suffix: $session")
        continue
    end
    println("Processing session: $session with Probe $prb")
    session_path = session_path * "\\"
    data = prep_session_data(session_path, prb, CHUNK)

    # Remove the processed nostril columns from the keypoint data
    for KP in (data.KP_R1, data.KP_R4, data.KP_R1_Cut, data.KP_R4_Cut)
        drop_columns!(KP, NOSTRIL_COLS)
    end

    println("Prefitting Encoder")
    β_eng, Σ_eng = prefit_engaged_encoder(
        data.KP_R1, data.KP_R4, data.PCA_R1, data.PCA_R4, data.FCs_R1, data.FCs_R4;
        lags=LAGS, start_time=START_TIME, pre=PRE, post=POST)

    println("Fitting Model")
    model, lls = fit_switching_model(
        data.KP_R1_Cut, data.KP_R4_Cut, data.PCA_R1_Cut, data.PCA_R4_Cut, β_eng, Σ_eng;
        lags=LAGS, start_time=START_TIME, max_iters=MAX_ITERS)

    println("Calculating average state path and encoding accuracy")
    inference = compute_state_inference(
        model, data.KP_R1, data.KP_R4, data.PCA_R1, data.PCA_R4, data.Tongue_R1, data.Tongue_R4, data.time_vec;
        lags=LAGS, start_time=START_TIME, three_sec=THREE_SEC, two_sec=TWO_SEC)

    accuracy = compute_encoding_accuracy(
        model, inference.X_R1_kernel, inference.X_R4_kernel,
        inference.Y_R1_trimmed, inference.Y_R4_trimmed,
        inference.R1_States, inference.R4_States, session; n_pcs=N_PCS)

    save_session_results(OUTPUT_DIR, session_save, session, accuracy, inference; n_pcs=N_PCS)
end
