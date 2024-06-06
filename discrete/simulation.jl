using Plots

include("state.jl")
include("walk.jl")

function show_stds(stds, qstds, steps)
    xs = collect(0:steps)
    plot(xs, stds, label="Classical Random Walk", linecolor=:blue, xlabel="Steps", ylabel="Standard Deviation", legend=:topleft, framestyle=:box, dpi=300)
    plot!(xs, qstds, label="Quantum Random Walk simulation", linecolor=:red)
    plot!([0, steps], [0, 0.54 * steps], label="Quantum Random Walk (≈ 0.54r)", linecolor=:black)
    savefig("outputs/discrete_$(steps).png")
end

function run_Z_simulation(; steps=50)
    qs0 = QState(1.0, 0, 0)
    qstds = walk(qs0; steps=steps)
    qstds

    s0 = State(1.0, 0)
    stds = walk(s0; steps=steps)
    stds

    show_stds(stds, qstds, steps)
end

run_Z_simulation(; steps=50)
run_Z_simulation(; steps=10)
