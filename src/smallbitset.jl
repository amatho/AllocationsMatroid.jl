mutable struct SmallBitSet{T<:Unsigned} <: AbstractSet{Int}
    bits::T

    SmallBitSet() = new{UInt64}(0)
    SmallBitSet{T}() where {T<:Unsigned} = new{T}(zero(T))
end

SmallBitSet{T}(itr) where {T<:Unsigned} = union!(SmallBitSet{T}(), itr)

function Base.empty!(s::SmallBitSet)
    s.bits = 0
    return s
end
Base.empty(::SmallBitSet{T}, ::Type{Int}=Int) where {T} = SmallBitSet{T}()
Base.emptymutable(::SmallBitSet{T}, ::Type{Int}=Int) where {T} = SmallBitSet{T}()

function Base.copy!(dest::SmallBitSet{T}, src::SmallBitSet) where {T}
    dest.bits = T(src.bits)
    return dest
end
Base.copy(s::SmallBitSet{T}) where {T} = copy!(SmallBitSet{T}(), s)
Base.copymutable(s::SmallBitSet) = copy(s)

function Base.length(s::SmallBitSet)
    return Base.count_ones(s.bits)
end

function Base.iterate(s::SmallBitSet, state=0)
    index = _next_index(s.bits, state)
    return isnothing(index) ? index : (index, index)
end

function Base.push!(s::SmallBitSet, x::Integer)
    _in_bounds(s.bits, x) || throw(DomainError(x, "value does not fit in a SmallBitSet of this size"))
    s.bits |= _from_index(x)
    return s
end

function Base.push!(s::SmallBitSet, items::Integer...)
    for x in items
        push!(s, x)
    end
    return s
end

function Base.delete!(s::SmallBitSet, x::Integer)
    _in_bounds(s.bits, x) && (s.bits = s.bits & ~_from_index(x))
    return s
end

function Base.union!(s1::SmallBitSet, s2::SmallBitSet)
    s1.bits |= s2.bits
    return s1
end

function Base.union!(s::SmallBitSet, itr)
    for x in itr
        push!(s, x)
    end
    return s
end

function Base.union!(s::SmallBitSet, r::AbstractUnitRange{<:Integer})
    a, b = first(r), last(r)
    a, b = max(a, 1), min(b, _typebits(s.bits))
    diff = b - a
    bits = (~zero(s.bits) >>> (_typebits(s.bits) - 1 - diff)) << (a - 1)
    s.bits |= bits
    return s
end

Base.union(s::SmallBitSet, sets...) = union!(copy(s), sets...)

function Base.intersect!(s1::SmallBitSet, s2::SmallBitSet)
    s1.bits &= s2.bits
    return s1
end

Base.intersect(s1::SmallBitSet, s2::SmallBitSet) = intersect!(copy(s1), s2)

function Base.setdiff!(s1::SmallBitSet, s2::SmallBitSet)
    s1.bits &= ~s2.bits
    return s1
end

function Base.symdiff!(s1::SmallBitSet, s2::SmallBitSet)
    s1.bits = xor(s1.bits, s2.bits)
    return s1
end

Base.issubset(a::SmallBitSet, b::SmallBitSet) = a.bits == a.bits & b.bits
Base.:⊊(a::SmallBitSet, b::SmallBitSet) = a.bits <= b.bits && a.bits != b.bits

@inline Base.in(x::Int, s::SmallBitSet) = _in_bounds(s.bits, x) ? _get_index(s.bits, x) : false
@inline Base.in(x::Integer, s::SmallBitSet) = _in_bounds(s.bits, x) ? _get_index(s.bits, x) : false

Base.first(s::SmallBitSet) = _next_index(s.bits, 0)
Base.last(s::SmallBitSet) = _prev_index(s.bits, _typebits(s.bits) - 1)
Base.minimum(s::SmallBitSet) = first(s)
Base.maximum(s::SmallBitSet) = last(s)
Base.extrema(s::SmallBitSet) = (first(s), last(s))

Base.:(==)(x::SmallBitSet, y::SmallBitSet) = x.bits == y.bits
Base.isequal(x::SmallBitSet, y::SmallBitSet) = isequal(x.bits, y.bits)

@inline _in_bounds(t::Unsigned, x::Integer) = x > 0 && x <= _typebits(t)
@inline _typebits(::UInt8) = 8
@inline _typebits(::UInt16) = 16
@inline _typebits(::UInt32) = 32
@inline _typebits(::UInt64) = 64
@inline _typebits(::UInt128) = 128
@inline _from_index(x::Integer) = unsigned(1 << (x - 1))
@inline _get_index(bits::T, idx::Integer) where {T<:Unsigned} = bits & _from_index(idx) != 0

function _next_index(bits::T, start::Integer) where {T<:Unsigned}
    start < 0 && (start = 0)
    start >= _typebits(bits) && return nothing
    i = start
    bits >>>= i
    while bits != 0
        if bits & 1 == 1
            return i + 1
        end
        i += 1
        bits >>>= 1
    end
    return nothing
end

function _prev_index(bits::T, start::Integer) where {T<:Unsigned}
    start < 0 && return nothing
    start >= _typebits(bits) && (start = _typebits(bits) - 1)
    i = start
    bits <<= (_typebits(bits) - 1 - i)
    while bits != 0
        if bits & (1 << (_typebits(bits) - 1)) > 0
            return i + 1
        end
        i -= 1
        bits <<= 1
    end
    return nothing
end

function bits_to_set(bits::Integer)
    s = SmallBitSet()
    s.bits = bits
    return s
end
