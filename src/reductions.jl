"""
    reduceinstance(V::Additive, C::Array{OrderedCategory,1}, agent::Int, removedbundle::Set{Int})

Reduce the instance given by the pair (V, C) to a new instance by giving the
supplied agent the supplied bundle. Returns an additive weight matrix, a set of
ordered categories for the new reduced instance and a function that turns an
allocation in the reduced instance into one for the original instance,
including giving the supplied agent the supplied bundle.
"""
function reduce_instance(V::Additive, C::Array{OrderedCategory,1}, agent::Int, removedbundle::Set{Int})
    N, M = agents(V), items(V)
    n, m = na(V), ni(V)

    Vs = zeros(n - 1, m - length(removedbundle))
    translate = zeros(Int, m - length(removedbundle))

    itemcounterchange = 0

    for j in M
        if j in removedbundle
            itemcounterchange += 1
            continue
        end

        newj = j - itemcounterchange
        translate[newj] = j

        Vs[1:n-1,newj] = [value(V, i, j) for i in N if i != agent]
    end

    # Create new ordered categories
    Cs = OrderedCategory[]
    index = 1
    for category in C
        newlength = length(category ∩ translate)
        push!(Cs, OrderedCategory(index, newlength, category.threshold))
        index += newlength
    end

    return Additive(Vs), Cs, (allocation) -> revert_instance(translate, agent, removedbundle, allocation)
end


"""
    revert_instance(translate::Array{Int, 1}, agent::Int, removedbundle::Set{Int}, allocation::Array{Set{Int}, 1})

Convert an allocation in a reduced instance to one in the original instance,
including giving the removed bundle to the removed agent.
"""
function revert_instance(translate::Array{Int, 1}, agent::Int, removedbundle::Set{Int}, allocation::Array{Set{Int}, 1})
    allocation = Set{Int}[Set{Int}(translate[item] for item in bundle) for bundle in allocation]
    insert!(allocation, agent, removedbundle)
end


"""
    create_ordered_instance(V::Additive, C::Counts)

Create an ordered instance for the given weights and categories. The items are
reorded such that each category has a continous range of indices for its items.
Returns a new additive weight matrix, an array of `OrderedCategory` objects and
a function that converts an allocation in the ordered instance to one in the
original instance.
"""
function create_ordered_instance(V::Additive, C::Counts)
    N = agents(V)
    Vo = zeros(na(V), ni(V))
    Co = OrderedCategory[]

    itemcounter = 1
    for category in C
         for i in N
            Vo[i,itemcounter:itemcounter+length(category)-1] = sort([value(V, i, j) for j in category], rev=true)
        end

        push!(Co, OrderedCategory(itemcounter, length(category), threshold(category)))
        itemcounter += length(category)
    end

    return Additive(Vo), Co, (alloc) -> revert_to_non_ordered_instance(V, C, Co, alloc)
end


"""
    revert_to_non_ordered_instance(V::Additive, C::Counts, Co::Array{OrderedCategory, 1}, alloc::Array{Set{Int}})

Convert an allocation in the ordered instance to one in the original instance.
"""
function revert_to_non_ordered_instance(V::Additive, C::Counts, Co::Array{OrderedCategory, 1}, alloc::Array{Set{Int}})
    bundles = [Set{Int}() for i in agents(V)]
    translate = zeros(Int, ni(V)) 

    # Create reverse lookup for items and agents
    for (i, bundle) in enumerate(alloc)
        for j in bundle
            translate[j] = i
        end
    end
 
    for (orig, new) in zip(C, Co)
        items = copy(orig.members)
        for j in new 
            i = translate[j]
            item = maximum((el) -> (value(V, i, el), el), items)
            push!(bundles[i], item[2])
            items = setdiff(items, item[2])
        end
    end

    return bundles
end

