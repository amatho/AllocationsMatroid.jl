"""
Tests for the axioms for the closed sets of a matroid, as given by Knuth (1974).
"""


"""
The ground set is closed. E ∈ F.
"""
function matroid_c1(m)
    last(m.F) == Set((ground_set(m),))
end


"""
The intersection of two closed sets is a closed set. 
If A, B ∈ F, then A ∩ B ∈ F.
"""
function matroid_c2(m)
    F = reduce(∪, m.F)
    for A ∈ F
        for B ∈ F
            if !(intersect(A, B) in F)
                return false
            end
        end
    end
    return true
end


"""
If A ∈ F and a, b ∈ E - A, then b is a member of all sets containing A ∪ {a} if and only if a is a member of all sets containing A ∪ {b}.
"""
function matroid_c3(m)
    E = ground_set(m)
    F = reduce(∪, m.F)
    delete!(F, BitSet())
    for A ∈ F
        t1 = setdiff(E, A)
        while !isempty(t1)
            a = minimum(t1)
            t2 = setdiff(t1, a)
            while !isempty(t2)
                b = minimum(t2)
                ā = reduce(intersect, [B for B in F if issubset(union(A, a), B)])
                b̄ = reduce(intersect, [B for B in F if issubset(union(A, b), B)])

                if issubset(b, ā) != issubset(a, b̄)
                    return false
                end

                delete!(t2, b)
            end
            delete!(t1, a)
        end
    end
    return true
end
