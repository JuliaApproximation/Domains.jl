# Definition of the abstract Domain type and its interface

# The type Domain{T} is defined in IntervalSets.jl

eltype(::Type{<:Domain{T}}) where {T} = T
prectype(::Type{<:Domain{T}}) where {T} = prectype(T)
numtype(::Type{<:Domain{T}}) where {T} = numtype(T)

convert_numtype(d::Domain{T}, ::Type{U}) where {T,U} = convert(Domain{to_numtype(T,U)}, d)
convert_prectype(d::Domain{T}, ::Type{U}) where {T,U} = convert(Domain{to_prectype(T,U)}, d)

Domain(d) = convert(Domain, d)

# Concrete types can implement similardomain(d, ::Type{T}) where {T}
# to support convert(Domain{T}, d) functionality.
convert(::Type{Domain{T}}, d::Domain{T}) where {T} = d
convert(::Type{Domain{T}}, d::Domain{S}) where {S,T} = similardomain(d, T)

"Can the domains be promoted without throwing an error?"
promotable_domains(d1, d2) = promotable_eltypes(eltype(d1),eltype(d2))
promotable_eltypes(::Type{S}, ::Type{T}) where {S,T} =
    isconcretetype(promote_type(S, T))
promotable_eltypes(::Type{S}, ::Type{T}) where {S<:AbstractVector,T<:AbstractVector} =
    promotable_eltypes(eltype(S), eltype(T))

"Promote the given domains to have a common element type."
promote_domains() = ()
promote_domains(domains...) = promote_domains(domains)
promote_domains(domains) = convert_eltype.(mapreduce(eltype, promote_type, domains), domains)

promote_domains(domains::AbstractSet{<:Domain{T}}) where {T} = domains
promote_domains(domains::AbstractSet{<:Domain}) = Set(promote_domains(collect(domains)))

convert_eltype(::Type{T}, d::Domain) where {T} = convert(Domain{T}, d)
convert_eltype(::Type{T}, d) where {T} = _convert_eltype(T, d, eltype(d))
_convert_eltype(::Type{T}, d, ::Type{T}) where {T} = d
_convert_eltype(::Type{T}, d, ::Type{S}) where {S,T} =
    error("Don't know how to convert the `eltype` of $(d).")
# Some standard cases
convert_eltype(::Type{T}, d::AbstractArray) where {T} = convert(AbstractArray{T}, d)
convert_eltype(::Type{T}, d::AbstractRange) where {T} = map(T, d)
convert_eltype(::Type{T}, d::Set) where {T} = convert(Set{T}, d)

promote(d1::Domain, d2::Domain) = promote_domains((d1,d2))
promote(d1::Domain, d2) = promote_domains((d1,d2))
promote(d1, d2::Domain) = promote_domains((d1,d2))

"A `EuclideanDomain` is any domain whose eltype is `<:StaticVector{N,T}`."
const EuclideanDomain{N,T} = Domain{<:StaticVector{N,T}}

"A `VectorDomain` is any domain whose eltype is `Vector{T}`."
const VectorDomain{T} = Domain{Vector{T}}

const AbstractVectorDomain{T} = Domain{<:AbstractVector{T}}

CompositeTypes.Display.displaysymbol(d::Domain) = 'D'

"What is the Euclidean dimension of the domain?"
dimension(::Domain{T}) where {T} = euclideandimension(T)

"Is the given combination of point and domain compatible?"
iscompatiblepair(x, d) = _iscompatiblepair(x, d, typeof(x), eltype(d))
_iscompatiblepair(x, d, ::Type{S}, ::Type{T}) where {S,T} =
    _iscompatiblepair(x, d, S, T, promote_type(S,T))
_iscompatiblepair(x, d, ::Type{S}, ::Type{T}, ::Type{U}) where {S,T,U} = true
_iscompatiblepair(x, d, ::Type{S}, ::Type{T}, ::Type{Any}) where {S,T} = false
_iscompatiblepair(x, d, ::Type{S}, ::Type{Any}, ::Type{Any}) where {S} = true

# Some generic cases where we can be sure:
iscompatiblepair(x::SVector{N}, ::EuclideanDomain{N}) where {N} = true
iscompatiblepair(x::SVector{N}, ::EuclideanDomain{M}) where {N,M} = false
iscompatiblepair(x::AbstractVector, ::EuclideanDomain{N}) where {N} = length(x)==N

compatible_or_false(x, domain) =
    iscompatiblepair(x, domain) ? true : (@warn "`in`: incompatible combination of point: $(typeof(x)) and domain eltype: $(eltype(domain)). Returning false."; false)

compatible_or_false(x::AbstractVector, domain::AbstractVectorDomain) =
    iscompatiblepair(x, domain) ? true : (@warn "`in`: incompatible combination of vector with length $(length(x)) and domain '$(domain)' with dimension $(dimension(domain)). Returning false."; false)


"Promote point and domain to compatible types."
promote_pair(x, d) = _promote_pair(x, d, promote_type(typeof(x),eltype(d)))
_promote_pair(x, d, ::Type{T}) where {T} = convert(T, x), convert(Domain{T}, d)
_promote_pair(x, d, ::Type{Any}) = x, d
# Some exceptions:
# - matching types: avoid promotion just in case it is expensive
promote_pair(x::T, d::Domain{T}) where {T} = x, d
# - `Any` domain
promote_pair(x, d::Domain{Any}) = x, d
# - tuples: these are typically composite domains and the elements may be promoted later on
promote_pair(x::Tuple, d::Domain{<:Tuple}) = x, d
# - abstract vectors: promotion may be expensive
promote_pair(x::AbstractVector, d::AbstractVectorDomain) = x, d
# - SVector: promotion is likely cheap
promote_pair(x::AbstractVector{S}, d::EuclideanDomain{N,T}) where {N,S,T} =
    _promote_pair(x, d, SVector{N,promote_type(S,T)})


# At the level of Domain we attempt to promote the arguments to compatible
# types, then we invoke indomain. Concrete subtypes should implement
# indomain. They may also implement `in` in order to accept more types.
#
# Note that if the type of x and the element type of d don't match, then
# both x and the domain may be promoted (using convert(Domain{T}, d)) syntax).
in(x, d::Domain) = compatible_or_false(x, d) && indomain(promote_pair(x, d)...)


"""
Return a suitable tolerance to use for verifying whether a point is close to
a domain. Typically, the tolerance is close to the precision limit of the numeric
type associated with the domain.
"""
default_tolerance(d::Domain) = default_tolerance(prectype(d))
default_tolerance(::Type{T}) where {T <: AbstractFloat} = 100eps(T)


"""
`approx_in(x, domain::Domain [, tolerance])`

Verify whether a point lies in the given domain with a certain tolerance.

The tolerance has to be positive. The meaning of the tolerance, in relation
to the possible distance of the point to the domain, is domain-dependent.
Usually, if the outcome is true, it means that the distance of the point to
the domain is smaller than a constant times the tolerance. That constant may
depend on the domain.

Up to inexact computations due to floating point numbers, it should also be
the case that `approx_in(x, d, 0) == in(x,d)`. This implies that `approx_in`
reflects whether a domain is open or closed.
"""
approx_in(x, d::Domain) = approx_in(x, d, default_tolerance(d))

function compatible_or_false(x, d, tol)
    tol >= 0 || error("Tolerance has to be positive in `approx_in`.")
    compatible_or_false(x, d)
end

approx_in(x, d::Domain, tol) =
    compatible_or_false(x, d, tol) && approx_indomain(promote_pair(x, d)..., tol)

# Fallback to `in`
approx_indomain(x, d::Domain, tol) = in(x, d)


isapprox(d1::Domain, d2::Domain; kwds...) = d1 == d2

isreal(d::Domain) = isreal(eltype(d))

infimum(d::Domain) = minimum(d)  # if the minimum exists, then it is also the infimum
supremum(d::Domain) = maximum(d)  # if the maximum exists, then it is also the supremum


"""
Return a bounding box of the given domain.

A bounding box is an interval, a hyperrectangle or the full space. It is such that
each point in the domain also lies in the bounding box.
"""
boundingbox(d) = FullSpace{eltype(d)}()

"Return the boundary of the given domain as a domain."
function boundary end
const ∂ = boundary

"""
Return the normal of the domain at the point `x`.

It is assumed that `x` is a point on the boundary of the domain.
"""
function normal end

"""
Return the tangents of the domain at the point `x`. The tangents form a
basis for the tangent plane, perpendicular to the normal direction at `x`.
"""
function tangents end

# "Lazy representation of the boundary of a domain."
# struct Boundary{T,D} <: Domain{T}
#     domain  :: D
# end
#
# domain(d::Boundary) = d.domain
#
# Boundary(domain) = Boundary{eltype(domain)}(domain)
# Boundary{T}(domain) where {T} = Boundary{T,typeof(domain)}(domain)
# Boundary{T}(domain::Domain{T}) where {T} = Boundary{T,typeof(domain)}(domain)
# Boundary{T}(domain::Domain{S}) where {S,T} = Boundary{T}(convert(Domain{T}, domain))
#
# indomain(x, d::Boundary) = in(x, boundary(domain(d)))
