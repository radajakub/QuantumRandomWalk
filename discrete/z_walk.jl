using Plots

struct QState
    coefficient::Complex
    coin::Int
    vertex::Int

    QState(coefficient, coin, vertex) = new(complex(coefficient), coin, vertex)
end

struct State
    probability::Real
    vertex::Int
end

split(qs::QState, c1::Complex, c2::Complex) = [QState(c1 * qs.coefficient, 0, qs.vertex), QState(c2 * qs.coefficient, 1, qs.vertex)]
increase(qs::QState) = QState(qs.coefficient, qs.coin, qs.vertex + 1)
decrease(qs::QState) = QState(qs.coefficient, qs.coin, qs.vertex - 1)

# implements Hadamard coin flip
function coin_flip(qs::QState)
    c1 = Complex(1 / sqrt(2))
    c2 = Complex(0)
    if qs.coin == 0
        c2 = c1
    elseif qs.coin == 1
        c2 = -c1
    else
        error("Invalid coin value")
    end
    return split(qs, c1, c2)
end

# implements the shift operator
function shift(qs::QState)
    if qs.coin == 0
        return increase(qs)
    elseif qs.coin == 1
        return decrease(qs)
    else
        error("Invalid coin value")
    end
end

# reduce the state
function reduce(states::Vector{QState})
    dict = Dict{Tuple{Integer,Integer},Complex}()
    for s in states
        key = (s.coin, s.vertex)
        if haskey(dict, key)
            dict[key] += s.coefficient
        else
            dict[key] = s.coefficient
        end
    end
    new_states = Vector{QState}()
    for (key, value) in dict
        if abs(value) > 1e-10
            push!(new_states, QState(value, key[1], key[2]))
        end
    end
    return new_states
end

function reduce(states::Vector{State})
    dict = Dict{Integer,Real}()
    for s in states
        key = s.vertex
        if haskey(dict, key)
            dict[key] += s.probability
        else
            dict[key] = s.probability
        end
    end
    new_states = Vector{State}()
    for (key, value) in dict
        if abs(value) > 1e-10
            push!(new_states, State(value, key))
        end
    end
    return new_states
end

function step(states::Vector{QState})
    new_states = Vector{QState}()
    for s in states
        ss = shift.(coin_flip(s))
        append!(new_states, ss)
    end
    return reduce(new_states)
end

function step(states::Vector{State})
    new_states = Vector{State}()
    for s in states
        push!(new_states, State(s.probability / 2, s.vertex - 1))
        push!(new_states, State(s.probability / 2, s.vertex + 1))
    end
    return reduce(new_states)
end

function position_std(states::Vector{QState})
    expectation = sum([real(s.coefficient^2) * s.vertex for s in states])
    expectation2 = sum([real(s.coefficient^2) * s.vertex^2 for s in states])
    return sqrt(expectation2 - expectation^2)
end

function position_std(states::Vector{State})
    expectation = sum([s.probability * s.vertex for s in states])
    expectation2 = sum([s.probability * s.vertex^2 for s in states])
    return sqrt(expectation2 - expectation^2)
end

function walk(s0; steps=10)
    states = [s0]
    stds = [position_std(states)]
    for _ in 1:steps
        states = step(states)
        sigma = position_std(states)
        push!(stds, sigma)
    end
    return stds
end

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
