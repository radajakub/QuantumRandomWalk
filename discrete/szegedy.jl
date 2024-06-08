# implementation based on https://arxiv.org/pdf/1611.02238

using LinearAlgebra
using QuantumWalk
using LightGraphs

struct BipartiteDoubleCover
    N::Int
    V::Vector{Int}
    E::Vector{Tuple{Int,Int}}
    adj::Vector{Vector{Int}}

    dim::Int

    function BipartiteDoubleCover(edges::Vector{Tuple{Int,Int}})
        N = maximum(maximum.(edges))
        V = collect(1:N)

        adj = [[] for _ in 1:N]

        for e in edges
            push!(adj[e[1]], e[2])
            push!(adj[e[2]], e[1])
        end

        dim = 2^ceil(Int, log2(N))

        return new(N, V, edges, adj, dim)
    end
end

function Base.show(io::IO, bg::BipartiteDoubleCover)
    println(io, "Bipartite Double Cover with $(bg.N) vertices")
    println(io, "Original Graph:")
    for e in bg.E
        println(io, "  $(e[1]) -> $(e[2])")
    end
    println(io, "Bipartite Double Cover adjacency:")
    for x in 1:bg.N
        println(io, "  $(x) -> $(bg.adj[x])")
    end
end

deg(bg::BipartiteDoubleCover, v::Int) = length(bg.adj[v])

function tensor(bg::BipartiteDoubleCover, ones::Vector{Int})
    t = zeros(Float64, bg.dim)
    t[ones] .= 1.0
    return t
end

function tensor(bg::BipartiteDoubleCover, x::Int, y::Int)
    return kron(tensor(bg, [x]), tensor(bg, [y]))
end

phi(bg::BipartiteDoubleCover, x::Int) = (1 / sqrt(deg(bg, x))) * sum(tensor(bg, x, y) for y in bg.adj[x])

psi(bg::BipartiteDoubleCover, y::Int) = (1 / sqrt(deg(bg, y))) * sum(tensor(bg, x, y) for x in bg.adj[y])

function Rs(bg::BipartiteDoubleCover)
    phis = [phi(bg, x) for x in bg.V]
    psis = [psi(bg, y) for y in bg.V]
    R1 = 2 * sum(p * p' for p in phis) - I
    R2 = 2 * sum(p * p' for p in psis) - I
    return (R1, R2)
end

function Qs(bg::BipartiteDoubleCover, t::Int)
    In = I(bg.dim)
    t_tensor = tensor(bg, [t])
    Q1 = kron(In - 2 * t_tensor * t_tensor', In)
    Q2 = kron(In, In - 2 * t_tensor * t_tensor')
    return Q1, Q2
end

function step(bg::BipartiteDoubleCover, s::Vector{Float64}, t::Int)
    Q1, Q2 = Qs(bg, t)
    R1, R2 = Rs(bg)
    return R2 * Q2 * R1 * Q1 * s
end

function measure(bg::BipartiteDoubleCover, s::Vector{Float64}, t::Int)
    m = tensor(bg, [t])
    M = kron(m * m', I(length(m)))
    return s' * M * s
end

init(bg::BipartiteDoubleCover) = (1 / sqrt(bg.N * (bg.N - 1))) * sum([tensor(bg, x, y) for x in bg.V for y in bg.V if x != y])

function run_and_measure(bg::BipartiteDoubleCover; t::Int=2, steps::Int=1)
    s = init(bg)

    for _ in 1:steps
        s = step(bg, s, t)
    end

    return measure(bg, s, t)
end

function find_half(bg::BipartiteDoubleCover; t::Int=2)
    p = 0.0
    step = 0
    while p < 0.5
        p = run_and_measure(bg, t=t, steps=step)
        println("t=$(step) p(v=$(t)) = $p")
        step += 1
    end
end

function verify(bg::BipartiteDoubleCover; t::Int, steps=1)
    g = SimpleGraph(bg.N)
    for e in bg.E
        add_edge!(g, e[1], e[2])
    end

    search = QWSearch(Szegedy(g), [t])

    mine = run_and_measure(bg, t=t, steps=steps)
    # set same initial state as for my implementation
    ref = execute_single_measured(search, init(bg), steps)[t]
    println("Mine: $mine")
    println("Ref:  $ref")
end

# construct a graph (for compatibility I chose the number of vertices to be 2^k for some k)
#        3 ---- 6 
#      / |
# 1 - 2  4 ---- 8
#      \ |
#        5 ---- 7
bg = BipartiteDoubleCover([(1, 2), (2, 3), (2, 5), (3, 4), (3, 6), (4, 8), (4, 5), (5, 7)])

# run my implementation on target vertex t for steps steps 
run_and_measure(bg; t=2, steps=2)

# run the algorithm until the target vertex t is visited with probability greater than 0.5
find_half(bg; t=3)

# compare my implementation with the reference implementation for the same target vertex t and number of steps
verify(bg; t=1, steps=1)

