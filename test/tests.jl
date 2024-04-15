import JuMP
using Allocations
using Graphs
using Random: seed!, randperm
using Test

# For utilities tests:
using Allocations: bipartite_matching, lex_optimize!

# For Counts test:
using Allocations: Category

# For Matroids test:
using Allocations: matroid_partition_knuth73, find_shortest_path, exchange_graph


function runtests(; slow_tests = true)

    seed!(4252170447285279131)

@testset "Types" begin

    @testset "Allocations" begin

        @testset "Basics" begin

            n, m = 3, 5

            A = Allocation(n, m)

            @test string(A) == "[{}, {}, {}]"

            @test na(A) == length(agents(A)) == n
            @test ni(A) == length(items(A)) == m

            @test isempty(bundle(A, 2))

            give!(A, 2, 4)
            give!(A, 2, 3)
            give!(A, 1, 2)

            @test string(A) == "[{2}, {3, 4}, {}]"

            give!(A, 3, [1, 5])

            @test string(A) == "[{2}, {3, 4}, {1, 5}]"

            @test 4 in bundle(A, 2)

            deny!(A, 2, 4)

            @test string(A) == "[{2}, {3}, {1, 5}]"
            @test summary(A) == "Allocation with 3 agents and 5 items, " *
                                "1 unallocated"

            A₂ = Allocation(A)
            @test A == A₂
            @test A !== A₂

            A₃ = copy(A)
            @test A == A₃
            @test A !== A₃

            deny!(A₂, 2, 3)
            @test A != A₂
            @test A₂ != A₃

            deny!(A₃, 2, 3)
            @test A₂ == A₃

            deny!(A, 2, 3)
            @test A == A₂

            @test isempty(bundle(A, 2))

            @test length(A) == length(collect(A)) == 3

            @test !isempty(A)
            @test isempty(Allocation())

            for g = 1:m
                give!(A, 1, g)
            end

            @test summary(A) == "Allocation with 3 agents and 5 items"

            A = Allocation(2, 4)
            give!(A, 2, 1)
            fill_even!(A)
            @test owner(A, 2) == 1
            @test length(bundle(A, 1)) == length(bundle(A, 2)) == 2

        end

    end

    @testset "Profiles" begin

        let

            X = [1 2 3; 3 2 1]

            V = Additive(X)
            V′ = Profile(X)

            @test V == V′

            i, g, h = 2, 3, 2

            @test value(V, i, g) == 1

            A = Allocation(5, 10)

            @test value(V, i, A) == 0

            give!(A, i, g)
            give!(A, i, h)

            @test value(V, i, A) == 3

        end

        let

            n, m = 5, 10

            V = Additive(n, m)

            i, g = 2, 3

            @test length(agents(V)) == na(V) == n
            @test length(items(V)) == ni(V) == m
            @test value(V, i, g) == 0

            value!(V, i, g, 4)

            @test value(V, i, g) == 4

        end

    end

    @testset "Counts" begin

        C = Counts(
            [1, 2, 3] => 2,
            [4, 5, 6] => 1
        )

        @test C isa Counts
        for c in C
            @test c isa Category
        end

        @test C[2].threshold == 1

    end

    @testset "Conflicts" begin

        C = Conflicts(star_graph(10))

        @test C isa Conflicts

        @test C.graph isa AbstractGraph

    end

    @testset "Inclusion/exclusion" begin

        A = Allocation()

        C₁ = Forbidden(A)
        C₂ = Permitted(A)
        C₃ = Required(A)

        @test C₁ isa Forbidden
        @test C₁.alloc === A
        @test C₂ isa Permitted
        @test C₂.alloc === A
        @test C₃ isa Required
        @test C₃.alloc === A

    end

end

@testset "Utilities" begin

    @testset "Matching" begin

        X = [0 1 1 0 0 0
             0 0 0 0 0 0
             1 0 0 1 0 0
             0 0 1 0 0 0
             0 0 1 1 0 0
             0 0 0 0 0 1]

        # M = bipartite_matching(Bool.(X))
        M = falses(size(X))
        for (i, g) in bipartite_matching(X)
            M[i, g] = true
        end

        @test all(sum(M, dims=1) .<= 1)
        @test all(sum(M, dims=2) .<= 1)
        @test sum(M) == 5

    end

    @testset "Lexicographic optimization" begin

        model = JuMP.Model(Allocations.conf.MIP_SOLVER)

        JuMP.@variable(model, x >= 0)
        JuMP.@variable(model, y >= 0)

        objectives = [(JuMP.MOI.MIN_SENSE, -x), (JuMP.MOI.MAX_SENSE, y)]

        JuMP.@constraint(model, x + y <= 10)
        JuMP.@constraint(model, x <= 3)

        constraints = lex_optimize!(model, objectives, ϵ=0)

        @test JuMP.termination_status(model) in Allocations.conf.MIP_SUCCESS

        @test JuMP.value(x) ≈ 3
        @test JuMP.value(y) ≈ 7

        # Cleanup
        for con in constraints
            JuMP.delete(model, con)
        end

    end

end

@testset "Basic checks" begin

    A = Allocation(2, 2)

    give!(A, 1, 1)

    @test !check_partition(A)
    @test !check_complete(A)

    give!(A, 2, 2)

    @test check_partition(A)
    @test check_complete(A)

    give!(A, 1, 2)

    @test !check_partition(A)
    @test check_complete(A)

end

@testset "EF checks" begin

    n, m = 2, 3

    V = Additive(ones(n, m))

    A = Allocation(n, m)

    @test check_ef(V, A)

    give!(A, 1, 1)

    @test !check_ef(V, A)
    @test check_ef1(V, A)

    give!(A, 1, 2)
    give!(A, 2, 3)

    @test check_ef1(V, A)
    @test check_efx(V, A)

    value!(V, 2, 1, 2)

    @test check_ef1(V, A)
    @test !check_efx(V, A)

end

@testset "Measures" begin

    V = Additive([1 3; 3 1])
    A = Allocation(2, 2)
    give!(A, 1, 2)
    give!(A, 2, 1)

    @test nash_welfare(V, A) ≈ 9
    @test utility(V, A) ≈ 6
    @test prop_alpha(V, A) ≈ 3 / (4 / 2)

end

@testset "MIPs" begin

    V₀ = Additive(rand(1:10, 3, 15))

    @testset "MNW" begin

        V = V₀

        let res = alloc_mnw(V)

            @test res.alloc isa Allocation
            @test check_ef1(V, res.alloc)
            @test res.mnw > 0

        end

        let res = alloc_mnw([1 2 3; 4 3 1])

            @test string(res.alloc) == "[{3}, {1, 2}]"
            @test res.mnw ≈ 3 * (4 + 3)

            @test res.model isa JuMP.Model

        end

    end

    @testset "MNW with constraints" for C in [
            Conflicts(path_graph(ni(V₀))),
            Counts(
                [1, 2, 3, 4]     => 3,
                [5, 6, 7]        => 2,
                [8, 9, 10]       => 2,
                [11, 12, 13, 14] => 3,
                [15]             => 1
            )
        ]

        V = V₀

        res = alloc_mnw(V)
        resc = alloc_mnw(V, C)

        @test check(V, resc.alloc, C)

        @test resc.alloc isa Allocation
        @test resc.mnw > 0

        # Adding constraint can't improve objective.
        @test resc.mnw <= res.mnw

        @test res.model isa JuMP.Model
        @test resc.model isa JuMP.Model

    end

    @testset "MIP with inclusion/exclusion" begin

        V = Additive([1 1 2; 1 2 1])

        # Somewhat arbitrarily using MNW

        A₀ = Allocation(V)
        give!(A₀, 1, 1)
        A = alloc_mnw(V, Required(A₀)).alloc
        @test 1 in bundle(A, 1)

        A₀ = Allocation(V)
        give!(A₀, 1, 3)
        A = alloc_mnw(V, Forbidden(A₀)).alloc
        @test !(3 in bundle(A, 1))

        A₀ = Allocation(V)
        give!(A₀, 1, [1, 2])
        give!(A₀, 2, [1, 2, 3])
        A = alloc_mnw(V, Permitted(A₀)).alloc
        @test !(3 in bundle(A, 1))

        # Test that MMS doesn't work, because of asymmetry:
        @test_throws MethodError alloc_mms(V, Required(A₀))
        @test_throws MethodError alloc_mms(V, Forbidden(A₀))
        @test_throws MethodError alloc_mms(V, Permitted(A₀))

    end

    @testset "MIP with multiple constraints" begin

        V = Profile([1 2 3; 3 1 1])
        F = Forbidden(Allocation(V, 1 => 3))
        R = Required(Allocation(V, 2 => 2))
        C = Constraints(F, R)

        @test string(alloc_mnw(V).alloc) == "[{2, 3}, {1}]"
        @test string(alloc_mnw(V, C).alloc) == "[{1}, {2, 3}]"

    end

    @testset "EF1 with conflicts" begin

        nᵢ = 6
        n, m = nᵢ, 2nᵢ

        V = Additive(rand(n, m))

        # Random graph whose components have at most n items, to ensure EF1 is
        # possible.
        G = blockdiag(SimpleGraph(nᵢ, nᵢ + 1), SimpleGraph(nᵢ, nᵢ + 1))
        C = Conflicts(G)

        res = alloc_ef1(V, C)

        @test res.model isa JuMP.Model

        @test check_ef1(V, res.alloc)

        # Check that it's usable without constraints, even though we're not
        # supplying a single-argument MIP implementation:
        @test check_ef1(V, alloc_ef1(V, nothing).alloc)

    end

    @testset "MNW+EF1 with conflicts" begin

        # Specific example (V and G) where MNW does not lead to EF1:

        V = Additive([10  5  8  1  2  9; 6  1  9  6  8  9])

        G = SimpleGraph(ni(V))
        add_edge!(G, 3, 4)
        add_edge!(G, 5, 6)

        C = Conflicts(G)

        res = alloc_mnw(V, C)

        @test !check_ef1(V, res.alloc)

        res1 = alloc_mnw_ef1(V, C)

        @test res.model isa JuMP.Model
        @test res1.alloc isa Allocation
        @test check_ef1(V, res1.alloc)

        # MNW did not yield EF1, so enforcing EF1 should reduce MNW:
        @test res1.mnw < res.mnw

    end

    @testset "EFX" begin

        V = V₀

        alloc_efx([1 1 1; 1 1 1]) # Bug fix: EFX but not EF

        res = alloc_efx(V)

        A = res.alloc

        @test res.model isa JuMP.Model
        @test A isa Allocation

        @test check_partition(A)
        @test check_complete(A)

        @test check_efx(V, A)

    end

    @testset "Maximin" begin

        V = V₀

        res = alloc_mm(V)

        A = res.alloc
        N = agents(V)

        @test res.model isa JuMP.Model
        @test A isa Allocation
        @test res.mm == minimum(value(V, i, A) for i in N)

    end

    @testset "Lexicographic" begin

        # Sanity check of lexicographic optimization, using internal functions.

        V = Additive([1 0; 1 1])

        ctx = Allocations.init_mip(V, Allocations.conf.MIP_SOLVER)

        A = ctx.alloc_var

        ctx.objectives = [
            (JuMP.MOI.MAX_SENSE, value(V, 1, A))
            (JuMP.MOI.MAX_SENSE, value(V, 2, A))
        ]

        ctx = Allocations.solve_mip(ctx)

        @test string(ctx.alloc) == "[{1}, {2}]"

    end

    @testset "Leximin" begin

        for (case, expected, values) in [

                ([2 1; 1 2]         , "[{1}, {2}]"       , [2, 2])
                ([1 3; 1 2]         , "[{2}, {1}]"       , [1, 3])
                ([4 3; 3 2]         , "[{2}, {1}]"       , [3, 3])
                ([2 2 2 2; 1 1 1 1] , nothing            , [2, 4])
                ([3 3 3 3; 1 1 1 1] , nothing            , [3, 3])
                ([2 1 2 1; 1 2 1 2] , "[{1, 3}, {2, 4}]" , [4, 4])
                ([3 4 2 1; 5 5 1 1] , "[{2, 3}, {1, 4}]" , [6, 6])

                # From Feige, Sapir & Tauber, "A Tight Negative Example for MMS
                # Fair Allocations", WINE'22:
                ([1 16 23 26 4 10 12 19 9
                  1 16 22 26 4  9 13 20 9
                  1 15 23 25 4 10 13 20 9] , nothing     , [39, 40, 43])

            ]

            V = Profile(case)
            N = agents(V)

            res = alloc_lmm(V)
            A = res.alloc

            @test res.model isa JuMP.Model
            @test A isa Allocation

            isnothing(expected) || @test string(A) == expected

            @test sort([value(V, i, A) for i in N]) == values

        end

    end

    slow_tests &&
    @testset "MMS" begin

        V = V₀

        let res = alloc_mms(V)

            @test res.alloc isa Allocation
            @test mms_alpha(V, res.alloc, res.mmss) ≈ res.alpha

        end

        let res = alloc_mms([3 1 2; 4 4 5])

            @test res.model isa JuMP.Model
            @test length(res.mms_models) == 2
            @test all(isa.(res.mms_models, JuMP.Model))

            @test res.alpha ≈ 1.0
            @test res.mmss ≈ [3.0, 5.0]

        end

        let res = alloc_mms([2 1; 1 2])

            @test res.alpha ≈ 2.0

        end

        let res = alloc_mms([2 1; 1 2], cutoff=true)

            @test res.alpha ≈ 1.0

        end

    end

    @testset "GGI" begin

        V = V₀

        let res = alloc_ggi(V)

            @test res.alloc isa Allocation
            @test res.model isa JuMP.Model

        end

        let res = alloc_ggi([1 1 3; 1 1 2])

            @test string(res.alloc) == "[{3}, {1, 2}]"

        end

    end

    slow_tests &&
    @testset "limits, $func" for func in [
        alloc_ef1, alloc_efx, alloc_mnw, alloc_mnw_ef1, alloc_mm, alloc_ggi,
        alloc_mms]
        # Could add `cutoff=true` for `alloc_mms`

        V = V₀

        # `alloc_ef1` and `alloc_mnw_ef1` require a constraint, so we'll supply
        # `nothing` to all:
        C = nothing

        N, M = agents(V), items(V)
        n, m = na(V), ni(V)

        for (min_bundle, max_bundle, min_owners, max_owners) in [
                (nothing, 10, 2, nothing)
                (3, 3, nothing, nothing)
                (rand(1:2, n), rand(2:3, n), nothing, rand(1:3, m))
            ]

            lob = something(min_bundle, 0)
            hib = something(max_bundle, m)
            loo = something(min_owners, 0)
            hio = something(max_owners, n)

            if func == alloc_mms && !(allequal(lob) && allequal(hib))
                @test_throws AssertionError func(V, C,
                      min_bundle=min_bundle,
                      max_bundle=max_bundle,
                      min_owners=min_owners,
                      max_owners=max_owners)
                continue
            end

            A = func(V, C,
                  min_bundle=min_bundle,
                  max_bundle=max_bundle,
                  min_owners=min_owners,
                  max_owners=max_owners).alloc

            broadcast(N, lob, hib) do i, lo, hi
                @test lo <= length(bundle(A, i)) <= hi
            end

            broadcast(M, loo, hio) do i, lo, hi
                @test lo <= length(owners(A, i)) <= hi
            end

        end
    end

end

@testset "Algorithms" begin

    V₀ = Additive(rand(1:10, 4, 12))

    @testset "Random" begin

        V = V₀

        res = alloc_rand(V)

        A = res.alloc

        @test A isa Allocation
        @test check_partition(A)

    end

    @testset "Random with conflicts" begin

        n, m = 10, 50
        nv, ne = m, 30

        V = Additive(zeros(n, m)) # Irrelevant

        # Random graph
        # (Could fail the Δ < n check, if we're unlucky)
        C = Conflicts(SimpleGraph(nv, ne))

        res = alloc_rand(V, C)
        A = res.alloc

        @test A isa Allocation
        # Implied by check_partition, but useful checkpoint:
        @test check_complete(A)
        @test check_partition(A)
        @test check(V, A, C)

    end

    slow_tests &&
    @testset "BKV18(1)" begin

        α = 1.061 # Approximation ratio, for geomean version of NW

        for _ = 1:10

            n = rand(2:4)
            m = 3n
            # Identical valuations
            V = Profile(repeat(rand(1:5, 1, m), n))

            res = alloc_bkv18_1(V)
            A = res.alloc

            @test A isa Allocation
            @test check_partition(A)
            @test α*nash_welfare(V, A)^(1/n) ≥ alloc_mnw(V).mnw^(1/n)
            @test check_efx(V, A)

        end

    end

    @testset "BKV18(2)" begin

        V = Profile([1 1 0 1 0
                     0 1 0 0 1
                     1 1 1 0 1])

        res = alloc_bkv18_2(V)
        A = res.alloc

        @test A isa Allocation
        @test check_partition(A)
        @test res.mnw == alloc_mnw(V).mnw

        V = Profile([1 0; 1 0; 0 1])

        res = alloc_bkv18_2(V)
        A = res.alloc
        @test res.mnw == nash_welfare(V, A) == 1
        @test nash_welfare(V, A, nonzero=false) == 0

        @test alloc_bkv18_2(Profile([0 0; 0 0])).mnw == 0

        V = Profile(ones(2, 10))
        for _ = 1:10
            @test alloc_bkv18_2(V).mnw == 25
        end

        for _ = 1:10
            V = Profile(rand(Bool, rand(2:5), rand(5:15)))
            @test alloc_bkv18_2(V).mnw == alloc_mnw(V).mnw
        end

        # Regression test
        V = Profile([0 0 1 0 1 1; 0 1 0 1 1 1; 0 1 0 1 0 0])
        for _ = 1:10
            alloc_bkv18_2(V) # Should not throw an exception
        end

        # Even distribution of unvalued items:
        V = Profile([1 1 1 0 0 0 0 0 0
                     0 0 0 0 0 0 0 0 0
                     0 0 0 0 0 0 0 0 0])
        A = alloc_bkv18_2(V).alloc
        for g = 4:9
            @test !owned(A, g)
        end

        res = alloc_bkv18_2(V, complete=true)
        @test res.mnw == 3
        A = res.alloc
        for i in agents(A)
            @test length(bundle(A, i)) == 3
        end

    end

    @testset "GHSS18(4)" begin

        # A modified version of the valuation function from the instance
        # where the best possible allocation has α = 3/4. The MMS with this
        # valuation function and m = 2n items, is 2.
        function v(m, B)
            if length(B) != 2
                return length(B)
            end

            g, g′ = minimum(B), maximum(B)
            if g == 1 && g′ == m || floor(g / 2) == floor(g′ / 2)
                return 2
            end

            return 3/2
        end

        # A function that creates a new set, with all items being incremented by
        # 1 (if m ∈ B, then 1 is in the new set). Setting the valuation of an
        # agent to `v(m, inc(m, B))` gives the agent an MMS of 2.
        inc(m, B) = Set(mod(g, 1:m) for g in B)

        n = 5
        m = 2n
        Vf = vcat([B -> v(m, B) for _ in 1:(n - 1)], [B -> v(m, inc(m, B))])
        V = Submodular(Vf, m)

        # Test when the MMS of each agent is known
        res = alloc_ghss18_4(V, repeat([2], n))

        @test !res.fail

        # The value each agent receives should be at least 1/3 * μᵢ = 2/3
        A = res.alloc
        for i in agents(V)
            @test value(V, i, bundle(A, i)) ≥ 2/3
        end

        # Test when the MMS of each agent is unknown
        res = alloc_ghss18_4b(V)

        # The value each agent receives should be at least 1/3 * μᵢ = 2/3
        for i in agents(V)
            @test value(V, i, bundle(A, i)) ≥ 2/3
        end
    end

    slow_tests &&
    @testset "MMS approximation with card. constr." begin

        # A default test set for all algorithms
        V₁ = V₀
        C₁ = Counts(
            [1, 3, 7, 9]          => 1,
            [4, 6, 8, 10, 11, 12] => 3,
            [2, 5]                => 2,
        )
        MMSs₁ = [mms(V₁, i, C₁).mms for i in agents(V₁)]

        # Check if alg finds an allocation that is complete, does not break the
        # cardinality constraints and gives each agent a bundle valued at no
        # less than `α * MMSs[i]`.
        function test_alg(alg, α, V, C, MMSs)

            A = alg(V, C).alloc

            @test A isa Allocation
            @test check_partition(A)

            # The allocation must not break the cardinality constraints
            for i in agents(V), c in C
                @test sum(owner(A, g) == i for g in c) ≤ c.threshold
            end

            for i in agents(V)
                @test value(V, i, bundle(A, i)) ≥ α * MMSs[i]
            end

        end

        @testset "1/3-MMS - BB18(3)" begin

            test_alg(alloc_bb18_3, 1/3, V₁, C₁, MMSs₁)

        end

        @testset "1/2-MMS - HH22(1)" begin

            test_alg(alloc_hh22_1, 1/2, V₁, C₁, MMSs₁)


            # The second half of the bag filling algorithm, where the floor_n(C)
            # highest-valued items are not worth 1/2 and thus ceil_n(C) must be
            # used instead, does not often get run when a random instance is
            # created. Thus, this tests the workings of that part
            let

                C = Counts(
                    [1, 5, 6]           => 3,
                    [2]                 => 1,
                    [3]                 => 1,
                    [4]                 => 3,
                    [7, 8, 9, 10]       => 5,
                )

                V = Additive([
                    0.2 0.4 0.4 0.4 0.1 0.1 0.1 0.1 0.1 0.1;
                    0.2 0.4 0.4 0.4 0.1 0.1 0.1 0.1 0.1 0.1
                ])

                test_alg(alloc_hh22_1, 1/2, V, C, [1, 1])

            end

        end
    end

    slow_tests &&
    @testset "MMS approximation" begin

        @testset "2/3-MMS - GMT18" begin

            function checkvalidallocation(A, V)
                # Test that we have received an allocation
                @test A isa Allocation

                # Test that all items are allocated properly
                for g in items(V)
                    @test owner(A, g) isa Int
                end

                # Test that each agent receives at least (2/3)-MMS
                for i in agents(V)
                    @test value(V, i, bundle(A, i)) ≥ 2/3 * mms(V, i).mms
                end
            end

            @testset "Base test" begin
                V = V₀

                A = alloc_gmt18(V)

                checkvalidallocation(A, V)
            end

            @testset "All bundles contain an item valued at ≥ 2/3" begin
                V = Additive([
                    0.45 0.27 0.10;
                    0.49 0.49 0.49
                ])

                A = alloc_gmt18(V)

                checkvalidallocation(A, V)

                # Test that each agent has an item valued at (2/3)-MMS or higher
                for i in agents(V)
                    @test any(value(V, i, g) ≥ 2/3 * mms(V, i).mms
                                for g in bundle(A, i))
                end
            end

            # Checks issue #1:
            # https://github.com/mlhetland/Allocations.jl/issues/1
            @testset "No item valued at ≥ 1/3" begin
                V = Additive([
                    0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25
                    0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25
                ])

                A = alloc_gmt18(V)

                checkvalidallocation(A, V)
            end
        end

    end

end

@testset "Matroids" begin
    @testset "GraphicMatroid properties" begin
        G = smallgraph(:karate)
        M = GraphicMatroid(G)

        @test rank(M) == 33
        @test rank(M, []) == 0
        @test is_indep(M, [1,2,3])
        @test is_indep(M, [1,2,16])
        @test is_indep(M, [1,2,17]) == false
        @test is_closed(M, 1:50) == false
        @test is_closed(M, 1:78) == false

        # situation encountered during manual testing:
        g = SimpleGraph{Int64}(64, [[9, 10, 11, 12, 13, 14, 15, 16], [9, 10, 11, 12, 13, 14], [9, 10, 11, 12, 13, 14], [9, 10, 11, 16], [9, 10, 11, 12, 13, 14, 15, 16], [9], [9, 10, 12, 14, 15, 16], [9, 10, 11], [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15], [1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 14, 15], [1, 2, 3, 4, 5, 8, 9, 10, 12, 13, 14, 15], [1, 2, 3, 5, 7, 9, 10, 11, 13, 15, 16], [1, 2, 3, 5, 9, 10, 11, 12, 16], [1, 2, 3, 5, 7, 9, 10, 11, 15, 16], [1, 5, 7, 9, 10, 11, 12, 14, 16], [1, 4, 5, 7, 12, 13, 14, 15]])
        m = GraphicMatroid(g)
        A = Set([64, 61, 55, 29, 52, 12, 37, 19, 4, 6, 13, 45])

        for e in [33,41,21]
            @test is_indep(m, A ∪ e)
        end
    end

    @testset "UniformMatroid properties" begin
        U = UniformMatroid(10, 6)
        F = FreeMatroid(10)
        Z = ZeroMatroid(10)

        @test is_indep(U, 0)
        @test is_indep(F, 0)
        @test is_indep(Z, 0)
        @test rank(U) == 6
        @test rank(F) == 10
        @test rank(Z) == 0

        for i in 1:5
            S = randperm(10)[1:i]
            @test is_indep(U, S) || S
            @test is_indep(F, S) || S
            @test is_indep(Z, S) == false || S

            @test is_circuit(U, S) == false || S
            @test is_circuit(F, S) == false || S
            @test is_circuit(Z, S) == (i == 1) || S

            @test rank(U, S) == i || S
            @test rank(F, S) == i || S
            @test rank(Z, S) == 0 || S

            @test closure(U, S) == S || S
            @test closure(F, S) == S || S
            @test closure(Z, S) == Set(1:10) || S
        end

        S = randperm(10)[1:6]
        @test is_indep(U, S) == true || S
        @test is_circuit(U, S) == false|| S
        @test rank(U, S) == 6 || S
        @test closure(U, S) == Set(1:10) || S

        S = randperm(10)[1:7]
        @test is_indep(U, S) == false || S
        @test is_circuit(U, S) || S
        @test rank(U, S) == 6 || S
        @test closure(U, S) == Set(1:10) || S

        
        for i in 8:10
            S = randperm(10)[1:i]
            @test is_indep(U, S) == false || S
            @test is_circuit(U, S) == false || S
            @test rank(U, S) == 6 || S
            @test closure(U, S) == Set(1:10) || S
        end
    end

    @testset "Exchange graphs and transfer paths" begin
        # Every agent likes every item.
        matroids = [FreeMatroid(5) for _ in 1:5]
        
        # Every agent has 1 item.
        A = Allocation(5,5)
        for i in 1:5
            give!(A, i, i)
        end

        # Every item should be exchangeable with every other item.
        @test exchange_graph(matroids, A) == complete_digraph(5)

        # Zachary's Karate Club.
        k = smallgraph(:karate)
        @test find_shortest_path(k, [1], [16]) |> length == 4
        @test find_shortest_path(k, [1,2,3], [1,18,4]) == [1]

        # A graph with no edges has no transfer paths.
        g = SimpleGraph(10, 0)
        @test find_shortest_path(g, [1,2,3], [4,5,6]) === nothing
        @test find_shortest_path(g, [1,2,3,4], [4,5]) == [4]


        # Three small matroids that require some transfers on the exchange
        # graph to get a partition of the ground set.

        # The ground set:       E = [1       2       3       4       5]
        g1 = SimpleGraph{Int64}(5, [[2, 3], [1, 3], [1, 2], [4],    [5]])
        g2 = SimpleGraph{Int64}(5, [[1],    [3, 4], [2, 4], [2, 3], [5]])
        g3 = SimpleGraph{Int64}(5, [[1],    [2],    [4, 5], [3, 5], [3, 4]])

        ms = [GraphicMatroid(g1), GraphicMatroid(g2), GraphicMatroid(g3)]

        (partition, junk) = matroid_partition_knuth73(ms)
        for (i, set) in enumerate(partition)
            @test is_indep(ms[i], set) || "set $set not indep in ms[$i]"
        end
    end
end

return nothing

end
