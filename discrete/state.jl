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

