
"The composition of several maps."
struct ComposedMap{T,MAPS} <: CompositeLazyMap{T}
    maps    ::  MAPS
end

ComposedMap(maps...) = ComposedMap{domaintype(maps[1])}(maps...)
ComposedMap{T}(maps...) where {T} = ComposedMap{T,typeof(maps)}(maps)

# TODO: make proper conversion
similarmap(m::ComposedMap, ::Type{T}) where {T} = ComposedMap{T}(m.maps...)

codomaintype(m::ComposedMap) = codomaintype(m.maps[end])

# Maps are applied in the order that they appear in m.maps
applymap(m::ComposedMap, x) = applymap_rec(x, m.maps...)
applymap_rec(x) = x
applymap_rec(x, map1, maps...) = applymap_rec(map1(x), maps...)

# The size of a composite map depends on the first and the last map to be applied
# We check whether they are scalar_to_vector, vector_to_vector, etcetera
mapsize(m::ComposedMap) = _composed_mapsize(m, m.maps[end], m.maps[1], mapsize(m.maps[end]), mapsize(m.maps[1]))
_composed_mapsize(m, m_end, m1, S_end::Tuple{Int,Int}, S1::Tuple{Int,Int}) = (S_end[1],S1[2])
_composed_mapsize(m, m_end, m1, S_end::Tuple{Int,Int}, S1::Tuple{Int}) =
    is_vector_to_scalar(m_end) ? () : (S_end[1],)
_composed_mapsize(m, m_end, m1, S_end::Tuple{Int,Int}, S1::Tuple{}) =
    is_vector_to_scalar(m_end) ? () : (S_end[1],)
_composed_mapsize(m, m_end, m1, S_end::Tuple{Int}, S1::Tuple{Int,Int}) = (S_end[1],S1[2])
_composed_mapsize(m, m_end, m1, S_end::Tuple{Int}, S1::Tuple{Int}) = (S_end[1],)
_composed_mapsize(m, m_end, m1, S_end::Tuple{Int}, S1::Tuple{}) = (S_end[1],)
_composed_mapsize(m, m_end, m1, S_end::Tuple{}, S1::Tuple{Int,Int}) = (1,S1[2])
_composed_mapsize(m, m_end, m1, S_end::Tuple{}, S1::Tuple{Int}) = ()
_composed_mapsize(m, m_end, m1, S_end::Tuple{}, S1::Tuple{}) = ()

function jacobian(m::ComposedMap, x)
    f, fd = backpropagate(x, reverse(components(m))...)
    fd
end
backpropagate(x, m1) = (m1(x), jacobian(m1, x))
function backpropagate(x, m2, ms...)
    f, fd = backpropagate(x, ms...)
    m2(f), jacobian(m2, f) * fd
end

for op in (:inverse, :leftinverse, :rightinverse)
    @eval $op(cmap::ComposedMap) = ComposedMap(reverse(map($op, components(cmap)))...)
end

inverse(m::ComposedMap, x) = inverse_rec(x, reverse(components(m))...)
inverse_rec(x) = x
inverse_rec(x, map1, maps...) = inverse_rec(inverse(map1, x), maps...)

leftinverse(m::ComposedMap, x) = leftinverse_rec(x, reverse(components(m))...)
leftinverse_rec(x) = x
leftinverse_rec(x, map1, maps...) = leftinverse_rec(leftinverse(map1, x), maps...)

rightinverse(m::ComposedMap, x) = rightinverse_rec(x, reverse(components(m))...)
rightinverse_rec(x) = x
rightinverse_rec(x, map1, maps...) = rightinverse_rec(rightinverse(map1, x), maps...)

composedmap() = ()
composedmap(m) = m
composedmap(m1, m2) = composedmap1(m1, m2)
composedmap1(m1, m2) = composedmap2(m1, m2)
composedmap2(m1, m2) = ComposedMap(m1, m2)

composedmap(m1, m2, maps...) = composedmap(composedmap(m1, m2), maps...)

composedmap(m1::ComposedMap, m2::ComposedMap) =
    ComposedMap(components(m1)..., components(m2)...)
composedmap1(m1::ComposedMap, m2) = ComposedMap(components(m1)..., m2)
composedmap2(m1, m2::ComposedMap) = ComposedMap(m1, components(m2)...)

# Arguments to ∘ should be reversed before passing on to mapcompose
(∘)(map1::AbstractMap, map2::AbstractMap) = composedmap(map2, map1)


==(m1::ComposedMap, m2::ComposedMap) =
    ncomponents(m1) == ncomponents(m2) && all(map(isequal, components(m1), components(m2)))
hash(m::ComposedMap, h::UInt) = hashrec("ComposedMap", collect(components(m)), h)

Display.combinationsymbol(m::ComposedMap) = Display.Symbol('∘')
Display.displaystencil(m::ComposedMap) =
    composite_displaystencil(m; reversecomponents=true)
show(io::IO, mime::MIME"text/plain", m::ComposedMap) = composite_show(io, mime, m)
show(io::IO, m::ComposedMap) = composite_show_compact(io, m)

## Lazy multiplication

"The lazy multiplication of one or more maps."
struct MulMap{T,MAPS} <: CompositeLazyMap{T}
    maps    ::  MAPS
end

MulMap(maps::Map{T}...) where {T} = MulMap{T}(maps...)
MulMap{T}(maps::Map{T}...) where {T} = MulMap{T,typeof(maps)}(maps)
MulMap{T}(maps...) where {T} = _mulmap(T, convert.(Map{T}, maps)...)
_mulmap(::Type{T}, maps...) where {T} = MulMap{T,typeof(maps)}(maps)

similarmap(m::MulMap, ::Type{T}) where {T} = MulMap{T}(m.maps...)

applymap(m::MulMap, x) = reduce(*, applymap.(components(m), Ref(x)))

multiply_map() = ()
multiply_map(m) = m
multiply_map(m1, m2) = multiply_map1(m1, m2)
multiply_map1(m1, m2) = multiply_map2(m1, m2)
multiply_map2(m1, m2) = MulMap(m1, m2)

multiply_map(m1, m2, maps...) = multiply_map(multiply_map(m1, m2), maps...)

multiply_map(m1::MulMap, m2::MulMap) =
    MulMap(components(m1)..., components(m2)...)
multiply_map1(m1::MulMap, m2) = MulMap(components(m1)..., m2)
multiply_map2(m1, m2::MulMap) = MulMap(m1, components(m2)...)

function Display.displaystencil(m::MulMap)
    A = Any[]
    list = components(m)
    push!(A, "x -> ")
    push!(A, Display.SymbolObject(list[1]))
    push!(A, "(x)")
    for i in 2:length(list)
        push!(A, " * ")
        push!(A, Display.SymbolObject(list[i]))
        push!(A, "(x)")
    end
    A
end
show(io::IO, mime::MIME"text/plain", m::MulMap) = composite_show(io, mime, m)


## Lazy sum

"The lazy sum of one or more maps."
struct SumMap{T,MAPS} <: CompositeLazyMap{T}
    maps    ::  MAPS
end

SumMap(maps::Map{T}...) where {T} = SumMap{T}(maps...)
SumMap{T}(maps::Map{T}...) where {T} = SumMap{T,typeof(maps)}(maps)
SumMap{T}(maps...) where {T} = _summap(T, convert.(Map{T}, maps)...)
_summap(::Type{T}, maps...) where {T} = SumMap{T,typeof(maps)}(maps)

similarmap(m::SumMap, ::Type{T}) where {T} = SumMap{T}(m.maps...)

applymap(m::SumMap, x) = reduce(+, applymap.(components(m), Ref(x)))

sum_map() = ()
sum_map(m) = m
sum_map(m1, m2) = sum_map1(m1, m2)
sum_map1(m1, m2) = sum_map2(m1, m2)
sum_map2(m1, m2) = SumMap(m1, m2)

sum_map(m1, m2, maps...) = sum_map(sum_map(m1, m2), maps...)

sum_map(m1::SumMap, m2::SumMap) =
    SumMap(components(m1)..., components(m2)...)
sum_map1(m1::SumMap, m2) = SumMap(components(m1)..., m2)
sum_map2(m1, m2::SumMap) = SumMap(m1, components(m2)...)

function Display.displaystencil(m::SumMap)
    A = Any[]
    list = components(m)
    push!(A, "x -> ")
    push!(A, Display.SymbolObject(list[1]))
    push!(A, "(x)")
    for i in 2:length(list)
        push!(A, " + ")
        push!(A, Display.SymbolObject(list[i]))
        push!(A, "(x)")
    end
    A
end
show(io::IO, mime::MIME"text/plain", m::SumMap) = composite_show(io, mime, m)


# Define the jacobian of a composite map
jacobian(m::ComposedMap) = composite_jacobian(reverse(components(m))...)
composite_jacobian(map1) = jacobian(map1)
composite_jacobian(map1, map2) = multiply_map(jacobian(map1) ∘ map2, jacobian(map2))
function composite_jacobian(map1, map2, maps...)
    rest = ComposedMap(reverse(maps)..., map2)
    f1 = jacobian(map1) ∘ rest
    f2 = composite_jacobian(map2, maps...)
    multiply_map(f1, f2)
end

jacobian(m::MulMap) = mul_jacobian(components(m)...)
mul_jacobian() = ()
mul_jacobian(map1) = jacobian(map1)
mul_jacobian(map1, map2) = sum_map(multiply_map(jacobian(map1), map2), multiply_map(map1, jacobian(map2)))
function mul_jacobian(map1, map2, maps...)
    rest = multiply_map(map2, maps...)
    mul_jacobian(map1, rest)
end
function jacobian(m::MulMap, x)
    z = map(t -> applymap(t,x), components(m))
    zd = map(t -> jacobian(t, x), components(m))
    sum(prod(z[1:i-1]) * zd[i] * prod(z[i+1:end]) for i in 1:ncomponents(m))
end


jacobian(m::SumMap) = sum_jacobian(components(m)...)
sum_jacobian() = ()
sum_jacobian(map1) = jacobian(map1)
sum_jacobian(maps...) = sum_map(map(jacobian, maps)...)

jacobian(m::SumMap, x) = sum(jacobian(mc, x) for mc in components(m))
