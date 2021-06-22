
###########################
# An n-dimensional simplex
###########################

"Supertype of an N-dimensional simplex."
abstract type Simplex{T,C} <: Domain{T} end

isclosedset(::Simplex{T,:closed}) where {T} = true
isclosedset(::Simplex{T,:open}) where {T} = false

isopenset(d::Simplex) = !isclosedset(d)

"The unit simplex is a polytope with the origin and all unit vectors as vertices."
abstract type UnitSimplex{T,C} <: Simplex{T,C} end

const ClosedUnitSimplex{T} = UnitSimplex{T,:closed}
const OpenUnitSimplex{T} = UnitSimplex{T,:open}

UnitSimplex(n::Int) = DynamicUnitSimplex(n)
UnitSimplex(::Val{N}) where {N} = EuclideanUnitSimplex{N}()

UnitSimplex{T}(n::Int) where {T <: StaticTypes} = StaticUnitSimplex{T}(n)
UnitSimplex{T}(::Val{N}) where {N,T} = StaticUnitSimplex{T}(Val(N))
UnitSimplex{T}() where {T <: StaticTypes} = StaticUnitSimplex{T}()
UnitSimplex{T}(n::Int) where {T} = DynamicUnitSimplex{T}(n)

UnitSimplex{T,C}(n::Int) where {T <: StaticTypes,C} = StaticUnitSimplex{T,C}(n)
UnitSimplex{T,C}(::Val{N}) where {N,T,C} = StaticUnitSimplex{T,C}(Val(N))
UnitSimplex{T,C}() where {T <: StaticTypes,C} = StaticUnitSimplex{T,C}()
UnitSimplex{T,C}(n::Int) where {T,C} = DynamicUnitSimplex{T,C}(n)

insimplex_closed(x) = mapreduce( t-> t >= 0, &, x) && norm(x,1) <= 1
insimplex_open(x) = mapreduce( t-> t > 0, &, x) && norm(x,1) < 1
insimplex_closed(x, tol) = mapreduce( t-> t >= -tol, &, x) && norm(x,1) <= 1+tol
insimplex_open(x, tol) = mapreduce( t-> t > -tol, &, x) && norm(x,1) < 1+tol

indomain(x, ::ClosedUnitSimplex) = insimplex_closed(x)
indomain(x, ::OpenUnitSimplex) = insimplex_open(x)
approx_indomain(x, ::ClosedUnitSimplex, tolerance) = insimplex_closed(x, tolerance)
approx_indomain(x, ::OpenUnitSimplex, tolerance) = insimplex_open(x, tolerance)

isempty(::UnitSimplex) = false

==(d1::UnitSimplex, d2::UnitSimplex) =
    isclosedset(d1)==isclosedset(d2) && dimension(d1)==dimension(d2)
hash(d::UnitSimplex, h::UInt) = hashrec("UnitSimplex", isclosedset(d), dimension(d), h)

boundingbox(d::UnitSimplex{T}) where {T} = UnitCube{T}(dimension(d))

boundary(d::UnitSimplex{T,:open}) where {T} = EmptySpace{T}()


distance_to(d::UnitSimplex, x) = x ∈ d ? zero(prectype(d)) : minimum(distance_to(el, x) for el in components(boundary(d)))

function normal(d::UnitSimplex, x)
    z = similar(x)
    fill!(z, 0)
    if sum(x) ≈ 1
        fill!(z, one(eltype(z))/sqrt(length(z)))
    else
        index = findmax(x)[2]
        z[index] = -1
    end
    return convert(eltype(d), z/norm(z))
end


"A unit simplex whose dimension is determined by its element type."
struct StaticUnitSimplex{T,C} <: UnitSimplex{T,C}
end

StaticUnitSimplex(::Val{N}) where {N} = StaticUnitSimplex{SVector{N,Float64}}()

StaticUnitSimplex{T}() where {T} = StaticUnitSimplex{T,:closed}()

StaticUnitSimplex{T}(n::Int) where {T} =
    (@assert n == euclideandimension(T); StaticUnitSimplex{T}())
StaticUnitSimplex{T}(::Val{N}) where {N,T} =
    (@assert N == euclideandimension(T); StaticUnitSimplex{T}())

StaticUnitSimplex{T,C}(n::Int) where {T,C} =
    (@assert n == euclideandimension(T); StaticUnitSimplex{T,C}())
StaticUnitSimplex{T,C}(::Val{N}) where {N,T,C} =
    (@assert N == euclideandimension(T); StaticUnitSimplex{T,C}())

const EuclideanUnitSimplex{N,T,C} = StaticUnitSimplex{SVector{N,T},C}

EuclideanUnitSimplex{N}() where {N} = EuclideanUnitSimplex{N,Float64}()

## A StaticUnitSimplex{<:Number} equals the interval [0,1]  (open or closed)
convert(::Type{Interval}, d::StaticUnitSimplex{T,:closed}) where {T <: Number} =
    UnitInterval{T}()
convert(::Type{Interval}, d::StaticUnitSimplex{T,:open}) where {T <: Number} =
    OpenInterval{T}(0, 1)

canonicaldomain(::Equal, d::StaticUnitSimplex{T}) where {T<:Number} = convert(Interval, d)

boundary(d::StaticUnitSimplex{T}) where {T<:Number} = boundary(convert(Interval, d))


"A unit simplex with vector elements with variable dimension determined by a field."
struct DynamicUnitSimplex{T,C} <: UnitSimplex{T,C}
    dimension   ::  Int

    DynamicUnitSimplex{T,C}(n::Int) where {T,C} = new(n)
    DynamicUnitSimplex{T,C}(n::Int) where {T<:StaticTypes,C} =
        (@assert n == euclideandimension(T); new(n))
end

DynamicUnitSimplex(n::Int) = DynamicUnitSimplex{Vector{Float64}}(n)
DynamicUnitSimplex{T}(n::Int) where {T} = DynamicUnitSimplex{T,:closed}(n)

dimension(d::DynamicUnitSimplex) = d.dimension

const VectorUnitSimplex{T,C} = DynamicUnitSimplex{Vector{T},C}

VectorUnitSimplex(dimension) = VectorUnitSimplex{Float64}(dimension)

show(io::IO, d::EuclideanUnitSimplex{N,Float64,:closed}) where {N} = print(io, "UnitSimplex(Val($(N)))")
show(io::IO, d::VectorUnitSimplex{Float64,:closed}) = print(io, "UnitSimplex($(dimension(d)))")


center(d::EuclideanUnitSimplex{N,T}) where {N,T} = ones(SVector{N,T})/N
center(d::VectorUnitSimplex{T}) where {T} = ones(T, dimension(d))/dimension(d)

corners(d::UnitSimplex) = vcat([origin(d)], [ unitvector(d, i) for i in 1:dimension(d)])

interior(d::EuclideanUnitSimplex{N,T}) where {N,T} = EuclideanUnitSimplex{N,T,:open}()
closure(d::EuclideanUnitSimplex{N,T}) where {N,T} = EuclideanUnitSimplex{N,T,:closed}()
interior(d::VectorUnitSimplex{T}) where {T} = VectorUnitSimplex{T,:open}(dimension(d))
closure(d::VectorUnitSimplex{T}) where {T} = VectorUnitSimplex{T,:closed}(dimension(d))


# We pick the center point, because it belongs to the domain regardless of
# whether it is open or closed.
point_in_domain(d::UnitSimplex) = center(d)

similardomain(d::StaticUnitSimplex{S,C}, ::Type{T}) where {S,T,C} =
    StaticUnitSimplex{T,C}()
similardomain(d::DynamicUnitSimplex{S,C}, ::Type{T}) where {S,T,C} =
    DynamicUnitSimplex{T,C}(d.dimension)


simplex_face_map(a::Number, b::Number, c::SVector{2}, d::SVector{2}) =
	AffineMap((d-c)/(b-a), c - (d-c)/(b-a)*a)
function simplex_face_map(a::Number, b::Number, c::Vector, d::Vector)
	@assert length(c) == length(d) == 2
	AffineMap((d-c)/(b-a), c - (d-c)/(b-a)*a)
end

function boundary(d::StaticUnitSimplex{SVector{2,T},:closed}) where {T}
    d0 = UnitInterval{T}()
	T0 = zero(T)
	T1 = one(T)
	maps = [
		simplex_face_map(T0, T1, SVector(T0,T0), SVector(T1,T0)),
		simplex_face_map(T0, T1, SVector(T1,T0), SVector(T0,T1)),
		simplex_face_map(T0, T1, SVector(T0,T1), SVector(T0,T0))
	]
	faces = map(m -> ParametricDomain(m, d0), maps)
	UnionDomain(faces)
end

# function boundary(d::StaticUnitSimplex{SVector{N,T},:closed}) where {N,T}
# 	left2 = infimum(d)
# 	right2 = supremum(d)
# 	d0 = UnitSimplex{SVector{N-1,T},:closed}()
# 	T0 = zero(T)
# 	T1 = one(T)
#
# 	map1 = cube_face_map(left1, right1, left2, right2, 1, left2[1])
# 	MAP = typeof(map1)
# 	maps = MAP[]
# 	for dim in 1:N
# 		push!(maps, cube_face_map(left1, right1, left2, right2, dim, left2[dim]))
# 		push!(maps, cube_face_map(left1, right1, left2, right2, dim, right2[dim]))
# 	end
# 	faces = map(m -> ParametricDomain(m, d_unit), maps)
# 	UnionDomain(faces)
# end
