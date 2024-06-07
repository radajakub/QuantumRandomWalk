using LinearAlgebra
using Plots
using LaTeXStrings

mutable struct KGraph
    N::Integer
    amplitudes::Matrix{Real}

    function KGraph(N::Integer)
        # initial amplitude of each vertex
        a = 1 / sqrt(N * (N - 1))

        # prepare the amplitude matrix
        amplitudes = fill(a, N, N)
        for i in 1:N
            amplitudes[i, i] = 0
        end

        new(N, amplitudes)
    end
end

Base.show(io::IO, G::KGraph) = show(io, G.amplitudes)

probs(G::KGraph) = sum(G.amplitudes .^ 2, dims=2)
probs(G::KGraph, v::Integer) = probs(G)[v]

function C!(g::KGraph, v::Integer)
    # (1) negate the outgoing amplitudes from target vertex v
    g.amplitudes[v, :] .*= -1

    # (2) swap the outgoing amplitudes with their opposites
    for i in 1:g.N
        if i != v
            g.amplitudes[i, v], g.amplitudes[v, i] = g.amplitudes[v, i], g.amplitudes[i, v]
        end
    end
end

function S!(g::KGraph)
    # (1) compute means of the graph amplitudes for every vertex
    mus = sum(g.amplitudes, dims=2) / (g.N - 1)
    # (2) reverse the outgoing amplitudes for each vertex by their mean
    for i in 1:g.N, j in 1:g.N
        if i != j
            # change compared to the slides because we want to reflect the value accross the mean so we must add the mean one more time
            g.amplitudes[i, j] = 2 * mus[i] - g.amplitudes[i, j]
        end
    end
end

function walk(N::Integer; v::Integer, steps::Integer=1000)
    g = KGraph(N)
    probs(g)
    ps = []
    for _ in 1:steps
        C!(g, v)
        S!(g)
        @assert isapprox(sum(probs(g)), 1)
        push!(ps, probs(g, v))
    end
    ps
end

function show_probabilities(ps, name, path)
    plt = plot(; title=name, xlabel="Number of steps", ylabel="Probability of being at a fixed vertex", framestyle=:box, dpi=300)
    xs = collect(1:length(ps))

    plot!(plt, xs, ps, linecolor=:blue, linewidth=2, legend=false)
    scatter!(plt, xs, ps, markersize=2, markercolor=:blue, legend=false)

    savefig(path)
end

v = 2
steps = 100

ps4 = walk(4, v=v, steps=steps)
show_probabilities(ps4, L"$K_{4}$", "outputs/k4.png")

ps10 = walk(10, v=v, steps=steps)
show_probabilities(ps10, L"$K_{10}$", "outputs/k10.png")

ps50 = walk(50, v=v, steps=steps)
show_probabilities(ps50, L"$K_{50}$", "outputs/k50.png")

ps100 = walk(100, v=v, steps=steps)
show_probabilities(ps100, L"$K_{100}$", "outputs/k100.png")

ps1000 = walk(1000, v=v, steps=steps)
show_probabilities(ps1000, L"$K_{1000}$", "outputs/k1000.png")
