using Distributions
using Graphs

# Note: Although `Uniform` represents a uniform distribution on [0,1], the
# implementation (https://github.com/JuliaStats/Distributions.jl/blob/a9b0e3c9)
# actually samples from [0,1), by using rand().


"""
    rand_additive(; n=2:10, m=n->2n:4n, v=(n, m)->Uniform(), rng=default_rng())
    rand_profile(; n=2:10, m=n->2n:4n, v=(n, m)->Uniform(), rng=default_rng())

Generate a random additive valuation profile (an `Additive` object) with the
number of agents and items, agents and values chosen using `rand` with the
given values as the first argument. Here `n` is used directly, while `m` should
be a *function of the number of agents*, which returns an argument for `rand`.
Similarly, `v` should be a function of both the number of agents and items.

The defaults for `n` and `m` are taken from [Hummel and
Hetland](https://arxiv.org/abs/2104.06280), who based them on based on
real-world data from [Caragiannis et al.](https://doi.org/10.1145/3355902), with
`n` at most 10, and the average `m/n` approximately 3.

The distribution for the values can be univariate, in which case it is used
independently for each value. For example, to generate Gaussian valuations, with
the [`Distributions`](https://github.com/JuliaStats/Distributions.jl) package,
use `v=(n, m)->Normal()`.

However, it can also be multivariate, in which case each sample should be a
vector of the same length as the number of items, representing the valuation
function of a single agent. For example, to get a `Dirichlet` distribution,
where an agent's values sum to `1`, you could use `v=(n, m)->Dirichlet(m, 2)`.

It is also possible to use a matrix-variate distribution, such as a matrix
normal distribution, where each sample should then be an `m` by `n` matrix, with
each column representing the valuation function of a single agent.

!!! warning

    Note that the matrix samples should be the *transpose* of the ones used in
    the resulting profile. This is to maintain consistency with the multivariate
    distributions, which produce column vectors.

`rand_profile` is an alias for `rand_additive`.
"""
function rand_additive(; n=2:10, m=n->2n:4n, v=(n, m)->Uniform(), rng=default_rng())
    nn = rand(rng, n)
    mm = rand(rng, m(nn))
    vv = v(nn, mm)
    X = rand(rng, nn, mm)
    rand!(rng, vv, X')
    return Additive(X)
end

const rand_profile = rand_additive


"""
    rand_conflicts_ws98(m; k=2:2:div(m, 2), β=Uniform(), rng=default_rng())
    rand_conflicts_ws98(V::Profile; ...)

Generate a random `Conflicts` contraint, whose underlying graph is constructed
according to the [Watts–Strogatz model](https://doi.org/10.1038/30918). The
keyword arguments `k` and `β` specify the possible values for the corresponding
parameters \$k\$ and \$\\beta\$, which are generated using `rand`. The defaults
are taken from [Hummel and Hetland](https://arxiv.org/abs/2104.06280). Note that
the parameter \$k\$ should be an even number, which Watts and Strogatz assume to
be much smaller than \$m\$.
"""
function rand_conflicts_ws98(m; k=2:2:div(m, 2), β=Uniform(), rng=default_rng())
    G = watts_strogatz(m, rand(rng, k), rand(rng, β), rng=rng)
    return Conflicts(G)
end

rand_conflicts_ws98(V::Profile; kwds...) = rand_conflicts_ws98(ni(V); kwds...)


"""
    rand_conflicts_er59(m, p=Uniform(), rng=default_rng())
    rand_conflicts_er59(V::Profile; ...)
    rand_conflicts(m; ...)
    rand_conflicts(m::Profile; ...)

Generate a random `Conflicts` contraint, whose underlying graph is constructed
according to the [Erdős–Rényi
model](https://doi.org/10.5486%2FPMD.1959.6.3-4.12). The keyword argument `p`
specifies the possible values for the corresponding parameter \$p\$, which is
generated using `rand`.

`rand_conflicts` is an alias for `rand_conflicts_er59`.
"""
function rand_conflicts_er59(m; p=Uniform(), rng=default_rng())
    G = erdos_renyi(m, rand(rng, p), rng=rng)
    return Conflicts(G)
end

rand_conflicts_er59(V::Profile; kwds...) = rand_conflicts_er59(ni(V); kwds...)
const rand_conflicts = rand_conflicts_er59


"""
    rand_conflicts_ba02(m; k=1:m, rng=default_rng())
    rand_conflicts_ba02(V::Profile; ...)

Generate a random `Conflicts` contraint, whose underlying graph is constructed
according to the [Barabási–Albert
model](https://arxiv.org/abs/cond-mat/0106096). The keyword argument `k`
specifies the possible values for the corresponding parameter \$k\$, which is
generated using `rand`.
"""
function rand_conflicts_ba02(m; k=1:m, rng=default_rng())
    G = barabasi_albert(m, rand(rng, k), rng=rng)
    return Conflicts(G)
end

rand_conflicts_ba02(V::Profile; kwds...) = rand_conflicts_ba02(ni(V); kwds...)


"""
    rand_matroid_er59(m; n=nothing, rng=default_rng())
    rand_matroid_er59(V::Profile; ...)

Generate a random `GraphicMatroid`, whose underlying graph is constructed
according to the [Erdős–Rényi model
](https://doi.org/10.5486%2FPMD.1959.6.3-4.12). The keyword argument `n`
specifies the possible values for the corresponding parameter \$n\$, which is
generated using `rand`.
"""
function rand_matroid_er59(m; n=nothing, rng=default_rng())
    if isnothing(n)
        # Find the minimum vertices in a graph with `m` edges
        min_verts = ceil(Int, sqrt(2 * m) + (1 / 2))
        max_verts = ceil(Int, sqrt(2) * m)
        nn = rand(rng, min_verts:max_verts)
    else
        nn = rand(rng, n)
    end
    G = erdos_renyi(nn, m, rng=rng)
    return GraphicMatroid(G)
end

rand_matroid_er59(V::Profile; kwds...) = rand_matroid_er59(ni(V), kwds...)


"""
    knuth_matroid(m, X)
    knuth_matroid(V::Profile, X)

Knuth's matroid construction (1974). Generates a matroid in terms of its closed
sets, given by the size of the universe `m` and a list of enlargements `X`.
"""
function knuth_matroid(m, X)
    r = 1 # r is current rank +1 due to 1-indexing.
    F = [Set([BitSet()])]
    E = BitSet(1:m)
    rank = Dict{BitSet,Integer}(BitSet() => 0)

    while E ∉ F[r]
        # Initialize F[r+1].
        push!(F, Set{BitSet}())

        # Setup add_set.
        add_callback = x -> rank[x] = r
        add_function = x -> add_set!(x, F, r, rank, add_callback)

        generate_covers!(F, r, E, add_function)

        # Perform enlargements.
        if r <= length(X)
            if X[r] !== nothing
                for x in X[r]
                    add_function(x)
                end
            end
        end

        r += 1
    end

    return ClosedSetsMatroid(m, r - 1, F, rank)
end

knuth_matroid(V::Profile, X) = knuth_matroid(ni(V), X)


"""
    rand_matroid_knu74_1(m, P; rng=default_rng())
    rand_matroid_knu74_1(V::Profile, P; ...)

Randomized version of Knuth's matroid construction. Random matroids are
generated via the method of random coarsening described in the 1974 paper.
Accepts the size of the universe `m`, and a list `P = [p₁, p₂, …]`, where `pᵢ`
denotes the number of coarsenings to apply at rank `i`.
"""
function rand_matroid_knu74_1(m, P; rng=default_rng())
    r = 1
    F = [Set([BitSet()])]
    E = BitSet(1:m)
    rank = Dict{BitSet,Integer}(BitSet() => 0)

    # Setup add_set.
    add_callback = x -> rank[x] = r
    add_function = x -> add_set!(x, F, r, rank, add_callback)

    while E ∉ F[r]
        r > m && @warn "Rank is larger than universe!" m r

        # Initialize F[r+1].
        push!(F, Set{BitSet}())

        generate_covers!(F, r, E, add_function)

        # Perform coarsening.
        if r <= length(P)
            coarsen!(F, r, E, P[r], add_function, rng=rng)
        end

        r += 1
    end

    return ClosedSetsMatroid(m, r - 1, F, rank)
end

rand_matroid_knu74_1(V::Profile, P; kwds...) = rand_matroid_knu74_1(ni(V), P, kwds...)


"""
    rand_matroid_knu74_2(m, P; rng=default_rng())
    rand_matroid_knu74_2(V::Profile, P; ...)

Randomized version of Knuth's matroid construction, that also keeps track of
independent sets. This entails keeping storing the whole power set of the ground
set in the rank table, and quickly becomes infeasible for values of `m` much
larger than 16.
"""
function rand_matroid_knu74_2(m, P; rng=default_rng())
    r = 1
    E = BitSet(1:m)
    rank = Dict{BitSet,Integer}(BitSet() => 0)

    F::Vector{Set{BitSet}} = [Set([BitSet()])]
    I::Vector{Set{BitSet}} = [Set([BitSet()])]

    while E ∉ F[r]
        r > m && @warn "Rank is larger than universe!" m r

        # Initialize F[r+1] and I[r+1].
        push!(F, Set{BitSet}())
        push!(I, Set{BitSet}())

        # Setup add_set.
        add_callback = x -> mark_independent_subsets!(x, I, r, length(x), rank)
        add_function = x -> add_set!(x, F, r, rank, add_callback)

        generate_covers!(F, r, E, add_function)

        # Perform coarsening.
        if r <= length(P)
            coarsen!(F, r, E, P[r], add_function, rng=rng)
        end

        r += 1
    end

    return FullMatroid(m, r - 1, F, I, Set(), rank)
end

rand_matroid_knu74_2(V::Profile, P; kwds...) = rand_matroid_knu74_2(ni(V), P, kwds...)


"""
    generate_covers!(F, r, E, insert_fn)

Generates minimal closed sets for rank r+1 and inserts them into F[r+1], using
the supplied insert_fn. This function should take one argument, the newly added
set.
"""
function generate_covers!(F, r, E, insert_fn)
    for y in F[r]
        t = setdiff(E, y)
        # Find all sets in F[r+1] that already contain y and remove excess
        # elements from t.
        for x in F[r+1]
            if issubset(y, x)
                setdiff!(t, x)
            end
            if isempty(t)
                break
            end
        end
        # Insert y ∪ a for all a ∈ t.
        while !isempty(t)
            a = minimum(t)
            insert_fn(union(y, a))
            setdiff!(t, a)
        end
    end
end


"""
    coarsen!(F, r, E, count, add_function; rng=default_rng())

Apply the specified number of coarsenings to the matroid.
"""
function coarsen!(F, r, E, count, add_function; rng=default_rng())
    for _ in 1:count
        if E ∈ F[r+1]
            return
        end
        A = rand(rng, F[r+1])
        t = setdiff(E, A)
        a = rand(rng, t)
        delete!(F[r+1], A)

        add_function(union(A, a))
    end
end


function add_set!(x, F, r, rank, callback)
    done = false

    while !done
        done = true

        for y in F[r+1]
            """
            When we are comparing a set X which we are investigating whether is a
            closed set, with a set Y which we know (based on the sets that have been
            added thus far) is a closed set, we can encounter three scenarios:

            1. X ∩ Y is a closed set (it is in our rank table) of rank < r
                - Move on. If this is the case for all Y ∈ F[r+1] we add X to F[r+1].

            2. X ∩ Y is a closed set (it is in our rank table) of rank == r
                - X and Y are not closed sets. Remove Y from F[r+1] and call this
                function on X ∪ Y.

            3. We have not seen X ∩ Y before (this happens when we have closed sets
                of lower rank but similar cardinality).

                a. If |X ∩ Y| < r, we know that the rank of X ∩ Y is < r. Move on.

                b. If |X ∩ Y| >= r, we need to see if X ∩ Y ⊆ Z for some Z of lower
                    rank.
                    If not, remove Y from F[r+1] and call this function on X ∪ Y.
            """

            x_and_y = intersect(x, y)
            if haskey(rank, x_and_y) && rank[x_and_y] < r
                continue
            elseif !haskey(rank, x_and_y)
                if length(x_and_y) < r
                    continue
                else
                    r´ = check_rank(x_and_y, r, F)
                    if !isnothing(r´)
                        rank[x_and_y] = r´
                        continue
                    end
                end
            end

            # x ∩ y has rank > r, replace with x ∪ y.
            delete!(F[r+1], y)
            x = union(x, y)
            done = false
            break
        end
    end

    push!(F[r+1], x)
    callback(x)
end


"""
    mark_independent_subsets!(x, I, r, c, rank)

Given a closed set x,
1. simply return if rank[x] < r (we've seen this already)
2. add it to I if |x| = r
3. recursively call this func on all x' ⊂ x st |x'| = |x| - 1
"""
function mark_independent_subsets!(x, I, r, c, rank)
    if haskey(rank, x) && rank[x] <= r
        return
    end
    if c == r
        push!(I[r+1], x)
    end
    rank[x] = r
    t = x
    while !isempty(t)
        v = setdiff(t, minimum(t))
        mark_independent_subsets!(setdiff(x, minimum(union(t, v))), I, r, c - 1, rank)
        t = v
    end
end


function check_rank(v, r, F)
    for (i, Fi) in enumerate(F[1:r])
        for z ∈ Fi
            if issubset(v, z)
                return i - 1
            end
        end
    end
    return nothing
end


"""
    knuth_matroid_erect(m, enlargements)

An improved version of KMC we are also finding the independent sets and circuits
of the matroid during generation.

This version assigns the Hamming weight of all subsets of E upfront. This is
infeasible for values of n much larger than 16.
"""
function knuth_matroid_erect(m, enlargements)
    # Initialize.
    r = 1
    mask = BitSet(1:m)
    mask_bits = 2^m - 1
    rank = Dict{BitSet,Int}()

    # Populate rank table with 100+cardinality for all subsets of E.
    rank[BitSet()] = 100
    k = 1
    while (k <= mask_bits)
        for i in 0:k-1
            rank[bits_to_set(k + i)] = rank[bits_to_set(i)] + 1
        end
        k = k + k
    end

    F = [Set([BitSet()])] # F[r] is the family of closed sets of rank r-1.
    I = [Set([BitSet()])] # I[r] is the family of independent sets of rank r-1.
    rank[BitSet()] = 0

    while mask ∉ F[r]
        push!(F, Set{BitSet}())
        push!(I, Set{BitSet}())

        # Generate minimal closed sets for rank r+1.
        for y in F[r] # y is a closed set of rank r.
            t = setdiff(mask, y) # The set of elements not in y.
            # Find all sets in F[r+1] that already contain y and remove excess elements from t.
            for x in F[r+1]
                if issubset(y, x)
                    setdiff!(t, x)
                end
            end
            # Insert y ∪ a for all a ∈ t.
            while !isempty(t)
                x = union(y, minimum(t))
                insert_set!(x, F, r, rank)
                setdiff!(t, x)
            end
        end

        # Enlarge (if any).
        if r <= length(enlargements) && enlargements[r] !== nothing
            for set in enlargements[r]
                insert_set!(set, F, r, rank)
            end
        end

        # Assign rank to sets and add independent ones to I.
        for M in F[r+1]
            mark!(M, I, r, rank)
        end

        # Next rank.
        r += 1
    end

    C = Set{BitSet}() # C is the set of circuits (minimal dependent sets) for M.
    k = 1
    while k <= mask_bits
        for i in 0:k-1
            ki_set = bits_to_set(k + i)
            i_set = bits_to_set(i)
            if rank[ki_set] == rank[i_set]
                push!(C, ki_set)
                unmark!(ki_set, rank[i_set] + 101, rank, mask)
            end
        end
        k += k
    end

    return FullMatroid(m, r - 1, F, I, C, rank)
end

knuth_matroid_erect(V::Profile, enlargements) = knuth_matroid_erect(ni(V), enlargements)


"""
    insert_set!(x, F, r, rank)

Inserts set x into F[r+1], but augments x if it is necessary to ensure no two
sets in F[r+1] have an intersection of rank greater than r.
"""
function insert_set!(x, F, r, rank)
    for y in F[r+1] # +1 since Julia is 1-indexed.
        if rank[intersect(x, y)] < r
            continue
        end

        # x ∩ y has rank > r, replace x and y with x ∪ y.
        delete!(F[r+1], y)
        insert_set!(union(x, y), F, r, rank)
        return
    end

    push!(F[r+1], x)
end


"""
    mark!(m, I, r, rank)

Given a closed set m, sets rank[m']=r for all subsets m' ⊆ m whose rank is not
already ≤ r, and adds m' to I if it is independent (that is, if its rank equals
its cardinality).
"""
function mark!(m, I, r, rank)
    if haskey(rank, m) && rank[m] <= r
        return
    end
    if rank[m] == 100 + r
        push!(I[r+1], m)
    end
    rank[m] = r
    t = m
    while !isempty(t)
        v = setdiff(t, minimum(t))
        mark!(setdiff(m, minimum(union(t, v))), I, r, rank)
        t = v
    end
end


function unmark!(m, card, rank, mask)
    if rank[m] < 100
        rank[m] = card
        t = setdiff(mask, m)
        while !isempty(t)
            v = setdiff(t, minimum(t))
            unmark!(union(m, minimum(union(t, v))), card + 1, rank, mask)
            t = v
        end
    end
end


function set_to_bits(S)
    if length(S) == 0
        return 0
    end

    reduce(+, (2^(x - 1) for x in S), init=0)
end


function bits_to_set(S::Integer)
    BitSet(i for (i, c) in enumerate(digits(S, base=2)) if c == 1)
end
