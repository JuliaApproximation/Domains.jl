
"A `ProductDomain` represents the cartesian product of other domains."
abstract type ProductDomain{T} <: CompositeDomain{T} end

composition(d::ProductDomain) = Product()

components(d::ProductDomain) = d.domains

==(d1::ProductDomain, d2::ProductDomain) = mapreduce(==, &, components(d1), components(d2))
hash(d::ProductDomain, h::UInt) = hashrec("ProductDomain", collect(components(d)), h)

isempty(d::ProductDomain) = any(isempty, components(d))
isclosedset(d::ProductDomain) = all(isclosedset, components(d))
isopenset(d::ProductDomain) = all(isopenset, components(d))

issubset(d1::ProductDomain, d2::ProductDomain) =
	compatibleproductdims(d1, d2) && all(map(issubset, components(d1), components(d2)))

volume(d::ProductDomain) = prod(map(volume, components(d)))

distance_to(d::ProductDomain, x) = sqrt(sum(distance_to(component(d, i), x[i])^2 for i in 1:ncomponents(d)))

compatibleproductdims(d1::ProductDomain, d2::ProductDomain) =
	dimension(d1) == dimension(d2) &&
		all(map(==, map(dimension, components(d1)), map(dimension, components(d2))))

Display.combinationsymbol(d::ProductDomain) = Display.Times()
Display.displaystencil(d::ProductDomain) = composite_displaystencil(d)
show(io::IO, mime::MIME"text/plain", d::ProductDomain) = composite_show(io, mime, d)
show(io::IO, d::ProductDomain) = composite_show_compact(io, d)

boundary_part(d::ProductDomain{T}, domains, i) where {T} =
	ProductDomain{T}(domains[1:i-1]..., boundary(domains[i]), domains[i+1:end]...)

boundary(d::ProductDomain) = productboundary(d)
productboundary(d) = productboundary(d, components(d))
productboundary(d, domains) =
	UnionDomain(boundary_part(d, domains, i) for i in 1:length(domains))
productboundary(d, domains::Tuple) =
	UnionDomain(tuple((boundary_part(d, domains, i) for i in 1:length(domains))...))

boundingbox(d::ProductDomain{T}) where {T} = ProductDomain{T}(map(boundingbox, components(d)))

infimum(d::ProductDomain) = toexternalpoint(d, map(infimum, components(d)))
supremum(d::ProductDomain) = toexternalpoint(d, map(supremum, components(d)))

interior(d::ProductDomain) = ProductDomain(map(interior, components(d)))
closure(d::ProductDomain) = ProductDomain(map(closure, components(d)))


VcatDomainElement = Union{Domain{<:Number},EuclideanDomain}

ProductDomain(domains...) = _ProductDomain(map(Domain, domains)...)
_ProductDomain(domains...) = TupleProductDomain(domains...)
_ProductDomain(domains::VcatDomainElement...) = VcatDomain(domains...)
ProductDomain(domains::AbstractVector) = VectorProductDomain(domains)
# To create a tuple product domain, invoke ProductDomain{T}. Here, we splat
# and this may end up creating a VcatDomain instead.
ProductDomain(domains::Tuple) = ProductDomain(domains...)

ProductDomain{T}(domains::Tuple) where {T} = ProductDomain{T}(domains...)
ProductDomain{T}(domains...) where {T} = _TypedProductDomain(T, domains...)
_TypedProductDomain(::Type{SVector{N,T}}, domains...) where {N,T} = VcatDomain{N,T}(domains...)
_TypedProductDomain(::Type{T}, domains...) where {T<:Vector} = VectorProductDomain{T}(domains...)
_TypedProductDomain(::Type{T}, domains...) where {T<:Tuple} = TupleProductDomain{T}(domains...)
_TypedProductDomain(::Type{T}, domains...) where {T} = TupleProductDomain{T}(domains...)

productdomain() = ()
productdomain(d) = d
productdomain(d1, d2, d3...) = productdomain(productdomain(d1, d2), d3...)

productdomain(d1, d2) = productdomain1(d1, d2)
productdomain1(d1, d2) = productdomain2(d1, d2)
productdomain2(d1, d2) = ProductDomain(d1, d2)

productdomain(d1::ProductDomain, d2::ProductDomain) =
	ProductDomain(components(d1)..., components(d2)...)
productdomain1(d1::ProductDomain, d2) = ProductDomain(components(d1)..., d2)
productdomain2(d1, d2::ProductDomain) = ProductDomain(d1, components(d2)...)

# Only override cross for variables of type Domain, it may have a different
# meaning for other variables (like the vector cross product)
cross(x::Domain...) = productdomain(x...)

^(d::Domain, n::Int) = productdomain(ntuple(i->d, n)...)

similardomain(d::ProductDomain, ::Type{T}) where {T} = ProductDomain{T}(components(d))

canonicaldomain(d::ProductDomain) = ProductDomain(map(canonicaldomain, components(d)))

mapto_canonical(d::ProductDomain) = ProductMap(map(mapto_canonical, components(d)))
mapfrom_canonical(d::ProductDomain) = ProductMap(map(mapfrom_canonical, components(d)))

for CTYPE in (Parameterization, Equal)
	@eval canonicaldomain(ctype::$CTYPE, d::ProductDomain) =
		ProductDomain(canonicaldomain.(Ref(ctype), components(d)))
	@eval mapto_canonical(ctype::$CTYPE, d::ProductDomain) =
		ProductMap(mapto_canonical.(Ref(ctype), components(d)))
	@eval mapfrom_canonical(ctype::$CTYPE, d::ProductDomain) =
		ProductMap(mapfrom_canonical.(Ref(ctype), components(d)))
end


"""
A `VcatDomain` concatenates the element types of its member domains in a single
static vector.
"""
struct VcatDomain{N,T,DIM,DD} <: ProductDomain{SVector{N,T}}
	domains	::	DD
end

VcatDomain(domains::Union{Vector,Tuple}) = VcatDomain(domains...)
function VcatDomain(domains...)
	T = numtype(domains...)
	N = sum(map(dimension, domains))
	VcatDomain{N,T}(domains...)
end

VcatDomain{N,T}(domains::Union{AbstractVector,Tuple}) where {N,T} = VcatDomain{N,T}(domains...)
function VcatDomain{N,T}(domains...) where {N,T}
	DIM = map(dimension,domains)
	VcatDomain{N,T,DIM}(convert_numtype.(domains, T)...)
end

VcatDomain{N,T,DIM}(domains...) where {N,T,DIM} =
	VcatDomain{N,T,DIM,typeof(domains)}(domains)

tointernalpoint(d::VcatDomain{N,T,DIM}, x) where {N,T,DIM} =
	convert_fromcartesian(x, Val{DIM}())
toexternalpoint(d::VcatDomain{N,T,DIM}, y) where {N,T,DIM} =
	convert_tocartesian(y, Val{DIM}())

"""
A `VectorProductDomain` is a product domain of arbitrary dimension where the
element type is a vector, and all member domains have the same element type.
"""
struct VectorProductDomain{V<:AbstractVector,DD<:AbstractVector} <: ProductDomain{V}
	domains	::	DD

	function VectorProductDomain{V,DD}(domains::DD) where {V,DD}
		@assert eltype(eltype(domains)) == eltype(V)
		new(domains)
	end
end

VectorProductDomain(domains::AbstractVector) =
	VectorProductDomain{Vector{eltype(eltype(domains))}}(domains)

VectorProductDomain{V}(domains::AbstractVector{<:Domain{T}}) where {T,V<:AbstractVector{T}} =
	VectorProductDomain{V,typeof(domains)}(domains)
function VectorProductDomain{V}(domains::AbstractVector) where {T,V<:AbstractVector{T}}
	Tdomains = convert.(Domain{T}, domains)
	VectorProductDomain{V}(Tdomains)
end

# Convenience: allow constructor to be called with multiple arguments, or with
# a container that is not a vector
VectorProductDomain(domains::Domain...) = VectorProductDomain(domains)
VectorProductDomain(domains) = VectorProductDomain(collect(domains))
VectorProductDomain{V}(domains::Domain...) where {V} = VectorProductDomain{V}(domains)
VectorProductDomain{V}(domains) where {V} = VectorProductDomain{V}(collect(domains))

# the dimension equals the number of composite elements
dimension(d::VectorProductDomain) = ncomponents(d)

tointernalpoint(d::VectorProductDomain, x) =
	(@assert length(x) == dimension(d); x)
toexternalpoint(d::VectorProductDomain, y) =
	(@assert length(y) == dimension(d); y)





"""
A `TupleProductDomain` is a product domain that concatenates the elements of
its member domains in a tuple.
"""
struct TupleProductDomain{T,DD} <: ProductDomain{T}
	domains	::	DD
end

TupleProductDomain(domains::Vector) = TupleProductDomain(domains...)
TupleProductDomain(domains::Domain...) = TupleProductDomain(domains)
TupleProductDomain(domains...) = TupleProductDomain(map(Domain, domains)...)
function TupleProductDomain(domains::Tuple)
	T = Tuple{map(eltype, domains)...}
	TupleProductDomain{T}(domains)
end

TupleProductDomain{T}(domains::Vector) where {T} = TupleProductDomain{T}(domains...)
TupleProductDomain{T}(domains...) where {T} = TupleProductDomain{T}(domains)
function TupleProductDomain{T}(domains::Tuple) where {T <: Tuple}
	Tdomains = map((t,d) -> convert(Domain{t},d), tuple(T.parameters...), domains)
	TupleProductDomain{T,typeof(Tdomains)}(Tdomains)
end
TupleProductDomain{T}(domains::Tuple) where {T} =
	TupleProductDomain{T,typeof(domains)}(domains)
