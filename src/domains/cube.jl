
"A `HyperRectangle` is the cartesian product of intervals."
abstract type HyperRectangle{T} <: ProductDomain{T} end

convert(::Type{HyperRectangle}, d::ProductDomain) =
	ProductDomain(map(t->convert(AbstractInterval, t), components(d)))

boundingbox(d::HyperRectangle) = d

"Compute all corners of the hyperrectangle."
function corners(d::HyperRectangle)
	left = infimum(d)
	right = supremum(d)
    N = length(left)
    corners = [similar(point_in_domain(d)) for i in 1:2^N]
    # All possible permutations of the corners
    for i=1:2^length(left)
        for j=1:N
            corners[i][j] = ((i>>(j-1))%2==0) ? left[j] : right[j]
        end
    end
    [convert(eltype(d), p) for p in corners]
end

# map the interval [a,b] to a cube defined by c (bottom-left) and d (top-right)
cube_face_map(a::Number, b::Number, c::SVector{2}, d::SVector{2}) = AffineMap((d-c)/(b-a), c - (d-c)/(b-a)*a)
function cube_face_map(a::Number, b::Number, c::Vector, d::Vector)
	@assert length(c) == length(d) == 2
	AffineMap((d-c)/(b-a), c - (d-c)/(b-a)*a)
end

# map the cube defined by a and b, to the cube defined by c and d,
# where the dimension of c and d is one larger, and one of the coordinates (dim) is fixed to dimval.
function cube_face_map(a::SVector{M}, b::SVector{M}, c::SVector{N}, d::SVector{N}, dim, dimval) where {N,M}
	@assert N == M+1
	T = promote_type(eltype(a),eltype(c))
	A = MMatrix{N,M,T}(undef)
	B = MVector{N,T}(undef)
	fill!(A, 0)
	fill!(B, 0)
	B[dim] = dimval
	for m = 1:dim-1
		# scalar map along the dimension "m"
		mapdim = interval_map(a[m], b[m], c[m], d[m])
		A[m,m] = mapdim.A
		B[m] = mapdim.b
	end
	for m = dim+1:N
		mapdim = interval_map(a[m-1], b[m-1], c[m], d[m])
		A[m,m-1] = mapdim.A
		B[m] = mapdim.b
	end
	AffineMap(SMatrix{N,M}(A), SVector{N}(B))
end

function cube_face_map(a::Vector, b::Vector, c::Vector, d::Vector, dim, dimval)
	M = length(a)
	N = length(c)
	@assert length(a)==length(b)
	@assert length(c)==length(d)
	@assert N == M+1
	T = promote_type(eltype(a),eltype(c))
	A = Matrix{T}(undef, N, M)
	B = Vector{T}(undef, N)
	fill!(A, 0)
	fill!(B, 0)
	B[dim] = dimval
	for m = 1:dim-1
		# scalar map along the dimension "m"
		mapdim = interval_map(a[m], b[m], c[m], d[m])
		A[m,m] = mapdim.A
		B[m] = mapdim.b
	end
	for m = dim+1:N
		mapdim = interval_map(a[m-1], b[m-1], c[m], d[m])
		A[m,m-1] = mapdim.A
		B[m] = mapdim.b
	end
	AffineMap(A, B)
end

# The boundary of a rectangle is a collection of mapped lower-dimensional rectangles.
# Dimension 2 is a special case, because the lower dimension is 1 where we use
# scalars instead of vectors
function boundary(d::HyperRectangle{SVector{2,T}}) where {T}
	left = infimum(d)
	right = supremum(d)
	x1 = left[1]; y1 = left[2]; x2 = right[1]; y2 = right[2]
	d_unit = UnitInterval{T}()
	maps = [
		cube_face_map(zero(T), one(T), SVector(x1,y1), SVector(x2,y1)),
		cube_face_map(zero(T), one(T), SVector(x2,y1), SVector(x2,y2)),
		cube_face_map(zero(T), one(T), SVector(x2,y2), SVector(x1,y2)),
		cube_face_map(zero(T), one(T), SVector(x1,y2), SVector(x1,y1))
	]
	faces = map(m -> ParametricDomain(m, d_unit), maps)
	UnionDomain(faces)
end

function boundary(d::HyperRectangle{SVector{N,T}}) where {N,T}
	left2 = infimum(d)
	right2 = supremum(d)
	d_unit = UnitCube{SVector{N-1,T}}()
	left1 = infimum(d_unit)
	right1 = supremum(d_unit)

	map1 = cube_face_map(left1, right1, left2, right2, 1, left2[1])
	MAP = typeof(map1)
	maps = MAP[]
	for dim in 1:N
		push!(maps, cube_face_map(left1, right1, left2, right2, dim, left2[dim]))
		push!(maps, cube_face_map(left1, right1, left2, right2, dim, right2[dim]))
	end
	faces = map(m -> ParametricDomain(m, d_unit), maps)
	UnionDomain(faces)
end

function boundary(d::HyperRectangle{Vector{T}}) where {T}
	if dimension(d) == 2
		left = infimum(d)
		right = supremum(d)
		x1 = left[1]; y1 = left[2]; x2 = right[1]; y2 = right[2]
		d_unit = UnitInterval{T}()
		maps = [
			cube_face_map(zero(T), one(T), [x1,y1], [x2,y1]),
			cube_face_map(zero(T), one(T), [x2,y1], [x2,y2]),
			cube_face_map(zero(T), one(T), [x2,y2], [x1,y2]),
			cube_face_map(zero(T), one(T), [x1,y2], [x1,y1])
		]
	else
		left2 = infimum(d)
		right2 = supremum(d)
		d_unit = UnitCube(dimension(d)-1)
		left1 = infimum(d_unit)
		right1 = supremum(d_unit)

		map1 = cube_face_map(left1, right1, left2, right2, 1, left2[1])
		MAP = typeof(map1)
		maps = MAP[]
		for dim in 1:dimension(d)
			push!(maps, cube_face_map(left1, right1, left2, right2, dim, left2[dim]))
			push!(maps, cube_face_map(left1, right1, left2, right2, dim, right2[dim]))
		end
	end
	faces = map(m -> ParametricDomain(m, d_unit), maps)
	UnionDomain(faces)
end


"A `Cube` is a hyperrectangle with equal side lengths in each dimension."
abstract type Cube{T} <: HyperRectangle{T} end

"A cube in a fixed N-dimensional Euclidean space."
const EuclideanCube{N,T} = Cube{SVector{N,T}}

"A cube with vector elements of variable length."
const VectorCube{T} = Cube{Vector{T}}

# This will cause a warning if vectors of different length are used with cubes
iscompatiblepair(x::AbstractVector, d::VectorCube) = length(x) == dimension(d)


"The unit cube is the domain `[0,1]^d`."
abstract type UnitCube{T} <: Cube{T} end

component(d::UnitCube, i::Int) =
    (1 <= i <= dimension(d) || throw(BoundsError); UnitInterval{numtype(d)}())

volume(d::UnitCube) = 1


"A unit cube that is specified by the element type `T`."
struct StaticUnitCube{T} <: UnitCube{T}
end

StaticUnitCube() = StaticUnitCube{SVector{3,Float64}}()
StaticUnitCube(::Val{N}) where {N} = StaticUnitCube{SVector{N,Float64}}()

StaticUnitCube{T}(n::Int) where {T} =
    (@assert n == euclideandimension(T); StaticUnitCube{T}())
StaticUnitCube{T}(::Val{N}) where {N,T} =
    (@assert N == euclideandimension(T); StaticUnitCube{T}())

similardomain(d::StaticUnitCube, ::Type{T}) where {T<:StaticTypes} =
    StaticUnitCube{T}()
similardomain(d::StaticUnitCube, ::Type{T}) where {T} =
    DynamicUnitCube{T}(dimension(d))

components(d::StaticUnitCube{SVector{N,T}}) where {N,T} =
    ntuple(x->UnitInterval{T}(), Val(N))

"The unit cube in a fixed N-dimensional space."
const EuclideanUnitCube{N,T} = StaticUnitCube{SVector{N,T}}

EuclideanUnitCube{N}() where {N} = EuclideanUnitCube{N,Float64}()

const UnitSquare{T} = EuclideanUnitCube{2,T}


"A unit cube whose dimension is specified by a field."
struct DynamicUnitCube{T} <: UnitCube{T}
    dimension   ::  Int

    DynamicUnitCube{T}(n::Int) where {T} = new(n)
    DynamicUnitCube{T}(n::Int) where {T<:StaticTypes} =
        (@assert n == euclideandimension(T); new(n))
end

DynamicUnitCube(n::Int) = DynamicUnitCube{Vector{Float64}}(n)

dimension(d::DynamicUnitCube) = d.dimension

components(d::DynamicUnitCube) = map(x->UnitInterval{numtype(d)}(), 1:dimension(d))

similardomain(d::DynamicUnitCube, ::Type{T}) where {T} =
    DynamicUnitCube{T}(d.dimension)
similardomain(d::DynamicUnitCube, ::Type{T}) where {T <: StaticTypes} =
    StaticUnitCube{T}()

"The unit cube with vector elements of a given dimension."
const VectorUnitCube{T} = DynamicUnitCube{Vector{T}}

VectorUnitCube(n::Int = 3) = VectorUnitCube{Float64}(n)
VectorUnitSquare() = VectorUnitCube(2)


UnitCube(n::Int) = DynamicUnitCube(n)
UnitCube(::Val{N} = Val(3)) where {N} = EuclideanUnitCube{N}()

UnitCube{T}(n::Int) where {T <: StaticTypes} = StaticUnitCube{T}(n)
UnitCube{T}(::Val{N}) where {N,T} = StaticUnitCube{T}(Val(N))
UnitCube{T}() where {T <: StaticTypes} = StaticUnitCube{T}()
UnitCube{T}(n::Int) where {T} = DynamicUnitCube{T}(n)

UnitCube(domains::UnitInterval...) = UnitCube(domains)
UnitCube(domains::NTuple{N,UnitInterval{T}}) where {N,T} =
    UnitCube{SVector{N,T}}(domains)
UnitCube(domains::SVector{N,UnitInterval{T}}) where {N,T} =
    UnitCube{SVector{N,T}}(domains)

UnitCube{T}(domains::UnitInterval...) where {T} = UnitCube{T}(domains)
UnitCube{T}(domain::NTuple{N,<:UnitInterval}) where {N,T<:SVector{N}} =
    StaticUnitCube{T}()
UnitCube{T}(domain::SVector{N,<:UnitInterval}) where {N,T<:SVector{N}} =
    StaticUnitCube{T}()

# Constructor: careful about ambiguities with FixedInterval arguments below
ProductDomain(domains::UnitInterval...) = UnitCube(domains...)
ProductDomain(domains::NTuple{N,<:UnitInterval}) where {N} = UnitCube(domains)
ProductDomain(domains::SVector{N,<:UnitInterval}) where {N} = UnitCube(domains)
ProductDomain(domains::AbstractVector{<:UnitInterval{T}}) where {T} =
    VectorUnitCube{T}(length(domains))
ProductDomain{T}(domains::UnitInterval...) where {N,S,T<:SVector{N,S}} =
	UnitCube{T}(domains...)
ProductDomain{T}(domains::NTuple{N,<:UnitInterval}) where {N,S,T<:SVector{N,S}} =
	UnitCube{T}(domains)
ProductDomain{T}(domains::SVector{N,<:UnitInterval}) where {N,S,T<:SVector{N,S}} =
	UnitCube{T}(domains)
ProductDomain{T}(domains::AbstractVector{<:UnitInterval}) where {S,T<:Vector{S}} =
	VectorUnitCube{S}(length(domains))

## Display:
show(io::IO, d::EuclideanUnitCube{3,Float64}) = print(io, "UnitCube()")
show(io::IO, d::EuclideanUnitCube{N,Float64}) where {N} = print(io, "UnitCube(Val($(N)))")
show(io::IO, d::UnitSquare{Float64}) = print(io, "UnitSquare()")
show(io::IO, d::UnitSquare{T}) where {T} = print(io, "UnitSquare{$(T)}()")
show(io::IO, d::VectorUnitCube{Float64}) = print(io, "UnitCube($(dimension(d)))")

# set the display stencils to [] to opt-out of composite display for the types above
Display.displaystencil(d::EuclideanUnitCube{N,Float64}) where {N} = []
Display.displaystencil(d::UnitSquare) = []
Display.displaystencil(d::VectorUnitCube{Float64}) = []
Display.object_parentheses(d::EuclideanUnitCube{N,Float64}) where {N} = false
Display.object_parentheses(d::UnitSquare) = false
Display.object_parentheses(d::VectorUnitCube{Float64}) = false



"An N-dimensional rectangle is the cartesian product of closed intervals."
struct Rectangle{T} <: HyperRectangle{T}
    a   ::  T
    b   ::  T

    function Rectangle{T}(a::S,b::S) where {S,T}
        @assert length(a)==length(b)
        new(a,b)
    end
end

dimension(d::Rectangle) = length(d.a)
component(d::Rectangle, i::Int) = ClosedInterval(d.a[i],d.b[i])
components(d::Rectangle) = map(ClosedInterval, d.a, d.b)

Rectangle(a, b) = Rectangle(promote(a,b)...)
Rectangle(a::T, b::T) where {T} = Rectangle{T}(a, b)
Rectangle(a::NTuple{N,T}, b::NTuple{N,T}) where {N,T} =
    Rectangle(SVector{N,T}(a), SVector{N,T}(b))

Rectangle(domains::Tuple) = Rectangle(domains...)
Rectangle(domains::ClosedInterval...) = Rectangle(promote_domains(domains)...)
Rectangle(domains::ClosedInterval{T}...) where {T} =
    Rectangle(map(infimum, domains), map(supremum, domains))
Rectangle(domains::AbstractVector{<:ClosedInterval}) =
    Rectangle(map(infimum, domains), map(supremum, domains))
Rectangle(domains::Domain...) =
    error("The Rectangle constructor expects two points or a list of intervals (closed).")

Rectangle{T}(domains::Tuple) where {T} = Rectangle{T}(domains...)
Rectangle{T}(domains::ClosedInterval...) where {T} =
	Rectangle{T}(map(infimum, domains), map(supremum, domains))
Rectangle{T}(domains::AbstractVector{<:ClosedInterval}) where {T} =
    Rectangle{T}(map(infimum, domains), map(supremum, domains))
Rectangle{T}(domains::Domain...) where {T} =
    error("The Rectangle constructor expects two points or a list of intervals (closed).")

ProductDomain(domains::ClosedInterval...) = Rectangle(domains...)
ProductDomain(domains::AbstractVector{<:ClosedInterval}) =
    Rectangle(map(infimum, domains), map(supremum, domains))
ProductDomain{T}(domains::ClosedInterval...) where {T} = Rectangle{T}(domains...)
ProductDomain{T}(domains::AbstractVector{<:ClosedInterval}) where {T} =
    Rectangle{T}(map(infimum, domains), map(supremum, domains))



"The N-fold cartesian product of a fixed interval."
struct FixedIntervalProduct{N,T,D} <: Cube{SVector{N,T}}
end

component(d::FixedIntervalProduct{N,T,D}, i::Int) where {N,T,D} =
    (1 <= i <= dimension(d) || throw(BoundsError); D())
components(d::FixedIntervalProduct{N,T,D}) where {N,T,D} =
    ntuple(x->D(), Val(N))

volume(d::FixedIntervalProduct{N,T,D}) where {N,T,D} = volume(D())^N

FixedIntervalProduct(domains::NTuple{N,D}) where {N,D <: FixedInterval} =
	FixedIntervalProduct{N,eltype(D),D}()
FixedIntervalProduct(domains::SVector{N,D}) where {N,D <: FixedInterval} =
	FixedIntervalProduct{N,eltype(D),D}()

const ChebyshevProductDomain{N,T} = FixedIntervalProduct{N,T,ChebyshevInterval{T}}
ChebyshevProductDomain(::Val{N}) where {N} = ChebyshevProductDomain{N}()
ChebyshevProductDomain{N}() where {N} = ChebyshevProductDomain{N,Float64}()

ProductDomain(domains::D...) where {D <: FixedInterval} =
	FixedIntervalProduct(domains)
ProductDomain(domains::NTuple{N,D}) where {N,D <: FixedInterval} =
	FixedIntervalProduct(domains)
ProductDomain(domains::SVector{N,<:FixedInterval}) where {N} =
	FixedIntervalProduct(domains)
ProductDomain{T}(domains::D...) where {N,S,T<:SVector{N,S},D <: FixedInterval} =
	FixedIntervalProduct(convert.(Ref(Domain{S}), domains))
ProductDomain{T}(domains::NTuple{N,D}) where {N,S,T<:SVector{N,S},D <: FixedInterval} =
	FixedIntervalProduct(convert.(Ref(Domain{S}), domains))
ProductDomain{T}(domains::SVector{N,<:FixedInterval}) where {N,S,T<:SVector{N}} =
	FixedIntervalProduct(domains)
