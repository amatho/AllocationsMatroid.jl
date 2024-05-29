abstract type Matroid end


struct ClosedSetsMatroid <: Matroid
    n::Integer # Size of universe
    r::Integer # Final rank (r == length(F)).
    F::Vector{Set{SmallBitSet}} # Closed sets by rank
end
Base.copy(M::ClosedSetsMatroid) = ClosedSetsMatroid(M.n, M.r, M.F)


struct FullMatroid <: Matroid
    n::Integer
    r::Integer
    F::Vector{Set{SmallBitSet}} # Closed sets by rank
    I::Vector{Set{SmallBitSet}} # Independent sets by rank
    C::Union{Nothing,Set{SmallBitSet}} # Circuits
end

Base.copy(M::FullMatroid) = FullMatroid(M.n, M.r, M.F, M.I, M.C)


struct GraphicMatroid <: Matroid
    g::Graph
    n::Integer
    r::Integer
end

GraphicMatroid(g::Graph) = GraphicMatroid(g, ne(g), length(kruskal_mst(g)))
Base.copy(M::GraphicMatroid) = GraphicMatroid(M.g, M.n, M.r)


struct UniformMatroid <: Matroid
    n::Integer
    r::Integer
end

Base.copy(M::UniformMatroid) = UniformMatroid(M.n, M.r)


struct PartitionMatroid <: Matroid
    n::Integer
    categories::Vector{UnitRange{Integer}}
end

Base.copy(M::PartitionMatroid) = PartitionMatroid(M.n, M.categories)


FreeMatroid(n) = UniformMatroid(n, n)
ZeroMatroid(n) = UniformMatroid(n, 0)


ground_set(M::ClosedSetsMatroid) = SmallBitSet{UInt64}(1:M.n)
ground_set(M::FullMatroid) = SmallBitSet{UInt64}(1:M.n)
ground_set(M::UniformMatroid) = BitSet(1:M.n)
ground_set(M::GraphicMatroid) = BitSet(1:M.n)


"""
    is_indep(M::Matroid, S)

Independence oracle. Determines whether S is independent in M.
"""
function is_indep end


"""
    is_indep(M::ClosedSetsMatroid, S::Integer)

Determines whether a given set S is independent in the matroid M, given by the
closed sets of M grouped by rank. Uses (I.1) in Greene (1989).
"""
function is_indep(M::ClosedSetsMatroid, S)
    t = length(S)

    if t == 0
        return true
    elseif t > M.r
        return false
    end

    for F in M.F[t]
        if issubset(S, F)
            return false
        end
    end

    return true
end


function is_indep(M::FullMatroid, S)
    card = length(S)
    if card == 0
        return true
    elseif card + 1 > length(M.I)
        return false
    end
    return S in M.I[card+1]
end


is_indep(M::UniformMatroid, S) = length(S) <= M.r


function is_indep(M::GraphicMatroid, S)
    edgelist = [e for (i, e) in enumerate(edges(M.g)) if i in S]
    subgraph, _vmap = induced_subgraph(M.g, edgelist)
    return !is_cyclic(subgraph)
end


"""
    rank(M::Matroid, S)

Returns the rank of the set S in M, ie. the size of the largest independent
subset of S.
"""
function rank end


rank(M::Matroid) = M.r


"""
    function rank(M::ClosedSetsMatroid, S)

Rank ``oracle''. Returns the rank of the set S in the matroid M.
"""
function rank(M::ClosedSetsMatroid, S)
    for (r, Fr) in enumerate(M.F), B ∈ Fr
        if issubset(S, B)
            return r - 1
        end
    end
end


function rank(M::FullMatroid, S)
    for (r, Fr) in enumerate(M.F), F ∈ Fr
        if issubset(S, F)
            return r - 1
        end
    end
end


rank(M::UniformMatroid, S) = min(length(S), M.r)


rank(M::GraphicMatroid, S) = length(S) > 0 ? length(minimal_spanning_subset(M, S)) : 0


"""
    is_circuit(M::Matroid, S)

Determines whether S is a circuit in M (ie. rank(M, S) = |S|-1).
"""
function is_circuit end


"""
    function is_circuit(M::ClosedSetsMatroid, S::Integer)

Determines whether a given set S is a circuit in the matroid M, given by the
closed sets of M. Uses (C.1) and (C.2) in Greene (1989).
"""
function is_circuit(M::ClosedSetsMatroid, S)
    t = length(S)

    for F in M.F[t] # (C.1) S ⊆ F for some F ∈ F_{t-1}.
        if issubset(S, F)
            @goto C2
        end
    end
    return false

    @label C2
    for F in M.F[t-1] # (C.2) |S ∩ F| ≤ r(F) for all F ∈ F_{t-2}.
        if length(intersect(S, F)) > t - 2
            return false
        end
    end

    return true
end


is_circuit(M::FullMatroid, S) = S in M.C


is_circuit(M::UniformMatroid, S) = M.r == length(S) - 1


function is_circuit(M::GraphicMatroid, S)
    return rank(M, S) == length(S) - 1
end


"""
    minimal_spanning_subset(M::Matroid, S)

Finds a minimal spanning subset of S in M. If S is the ground set of M, this
produces a basis of M.
"""
function minimal_spanning_subset end


"""
    function minimal_spanning_subset(M::ClosedSetsMatroid, A::Integer)

Algorithm 3.1 from Greene (1989). Given a matroid M = (E, F) and some subset A
of E, finds a minimal spanning subset of A. If A = E, this finds a basis of M.
If A is a basis, this finds A.
"""
minimal_spanning_subset(M::ClosedSetsMatroid, A) = _mss(M, 0, A)


minimal_spanning_subset(M::FullMatroid, A) = _mss(M, 0, A)


function _mss(M::Union{ClosedSetsMatroid,FullMatroid}, j::Integer, A)
    B = [intersect(F, A) for F in M.F[j+1] if length(intersect(F, A)) > j]

    while length(B) == 0
        if j >= length(A) - 1
            return A
        end

        j += 1
        B = [intersect(F, A) for F in M.F[j+1] if length(intersect(F, A)) > j]
    end

    _mss(M, j, setdiff(A, rand(reduce(union, B))))
end


minimal_spanning_subset(M::UniformMatroid, S) = throw("unimplemented")


"""
    minimal_spanning_subset(M::GraphicMatroid, S)

Uses Kruskal's algorithm to find a minimal spanning tree over M.G.
"""
function minimal_spanning_subset(M::GraphicMatroid, S)
    edgelist = [e for (i, e) in enumerate(edges(M.g)) if i in S]
    subgraph, _vmap = induced_subgraph(M.g, edgelist)
    return kruskal_mst(subgraph)
end


"""
    function minimal_spanning_subsets(M::ClosedSetsMatroid, A::Integer)

A modification of Algorithm 3.1 from Greene (1989) that finds all minimal
spanning subsets of A ⊆ E, given a matroid M = (E, F). If A = E, this finds the
bases of M.
"""
minimal_spanning_subsets(M::ClosedSetsMatroid, A) = _mss_all(M, 0, A)


minimal_spanning_subsets(M::FullMatroid, A) = _mss_all(M, 0, A)


function _mss_all(M, j::Integer, A)
    B = [intersect(F, A) for F in M.F[j+1] if length(intersect(F, A)) > j]

    while length(B) == 0
        if j >= length(A) - 1
            # Make sure to add A to a new set, without flattening (since A is
            # already a set)
            S = Set{AbstractSet{Int}}()
            push!(S, A)
            return S
        end

        j += 1
        B = [intersect(F, A) for F in M.F[j+1] if length(intersect(F, A)) > j]
    end

    bases = Set{AbstractSet{Int}}()
    T = reduce(union, B)
    while !isempty(T)
        x = minimum(T)
        union!(bases, _mss_all(M, j, setdiff(A, x)))
        setdiff!(T, x)
    end

    return bases
end


minimal_spanning_subsets(M::UniformMatroid, A) = throw("unimplemented")
minimal_spanning_subsets(M::GraphicMatroid, A) = throw("unimplemented")


"""
    bases(M::Matroid)

Finds the set of bases of M.
"""
function bases end


bases(M::ClosedSetsMatroid) = _mss_all(M, 0, ground_set(M))
bases(M::FullMatroid) = _mss_all(M, 0, ground_set(M))


bases(M::UniformMatroid) = throw("unimplemented")
bases(M::GraphicMatroid) = throw("unimplemented")


"""
    closure(M::Matroid, S)

Finds the closure of a set S in M, that, when given a set S ⊆ E, returns the set
of elements in x ∈ E such that x can be added to S with no increase in rank. It
returns the closed set of the same rank as S, that contains S.
"""
function closure end


function closure(M::ClosedSetsMatroid, S)
    for Fr in M.F, B in Fr
        if issubset(S, B)
            return B
        end
    end
end


function closure(M::FullMatroid, S)
    for Fr in M.F, B in Fr
        if issubset(S, B)
            return B
        end
    end
end


closure(M::UniformMatroid, S) = length(S) < M.r ? S : ground_set(M)


function closure(M::GraphicMatroid, S)
    edgelist = [e for (i, e) in enumerate(edges(M.g)) if i in S]
    _sg, vmap = induced_subgraph(M.g, edgelist)
    return [e for e in edges(M.g) if [e.src, e.dst] ⊆ vmap || e.src == e.dst]
end


is_closed(M::Matroid, S) = closure(M, S) == S
