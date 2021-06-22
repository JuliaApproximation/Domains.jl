
"""
A `MappedDomain` represents the mapping of a domain.

The map of a domain `d` under the mapping `y=f(x)` consists of all points `f(x)`
with `x ∈ d`. The characteristic function of a mapped domain is defined in
terms of the inverse map `g = inverse(f)`:
```
x ∈ m ⟺ g(x) ∈ d
```
"""
abstract type AbstractMappedDomain{T} <: SimpleLazyDomain{T} end

superdomain(d::AbstractMappedDomain) = d.domain

const MappedVectorDomain{T} = AbstractMappedDomain{Vector{T}}

tointernalpoint(d::AbstractMappedDomain, x) = inverse_map(d, x)
toexternalpoint(d::AbstractMappedDomain, y) = forward_map(d, y)

canonicaldomain(d::AbstractMappedDomain) = canonicaldomain(superdomain(d))
mapfrom_canonical(d::AbstractMappedDomain) = forward_map(d) ∘ mapfrom_canonical(superdomain(d))
mapto_canonical(d::AbstractMappedDomain) = mapto_canonical(superdomain(d)) ∘ inverse_map(d)


# TODO: check whether the map alters the dimension
dimension(d::MappedVectorDomain) = dimension(superdomain(d))

# TODO: check whether the map affects these properties
isempty(d::AbstractMappedDomain) = isempty(superdomain(d))
isopenset(d::AbstractMappedDomain) = isopenset(superdomain(d))
isclosedset(d::AbstractMappedDomain) = isclosedset(superdomain(d))

corners(d::AbstractMappedDomain) = [forward_map(d, x) for x in corners(superdomain(d))]

## I/O functionality

show(io::IO, mime::MIME"text/plain", d::AbstractMappedDomain) = composite_show(io, mime, d)
Display.displaystencil(d::AbstractMappedDomain) =
    map_stencil_broadcast(forward_map(d), superdomain(d))
Display.object_parentheses(d::AbstractMappedDomain) =
    Display.object_parentheses(forward_map(d))
Display.stencil_parentheses(d::AbstractMappedDomain) =
    Display.stencil_parentheses(forward_map(d))



"A `MappedDomain` stores the inverse map of a mapped domain."
struct MappedDomain{T,F,D} <: AbstractMappedDomain{T}
    invmap  ::  F
    domain  ::  D
end

# In the constructor, we have to decide which T to use for the MappedDomain.
# - we don't know anything about invmap: deduce T from the given domain
MappedDomain(invmap, domain::Domain{T}) where {T} = MappedDomain{T}(invmap, domain)
# - if the map is a Map{T}, use that T for the MappedDomain
MappedDomain(invmap::Map{T}, domain::Domain{S}) where {S,T} = MappedDomain{T}(invmap, domain)
# - if T is given in the constructor, by all means we use that
MappedDomain{T}(invmap, domain::Domain) where {T} =
    MappedDomain{T,typeof(invmap),typeof(domain)}(invmap, domain)
# - in that case, if the map is a Map{S}, make sure that S matches T
MappedDomain{T}(invmap::Map{T}, domain::Domain) where {T} =
    MappedDomain{T,typeof(invmap),typeof(domain)}(invmap, domain)
MappedDomain{T}(invmap::Map{S}, domain::Domain) where {S,T} =
    MappedDomain{T}(convert(Map{T}, invmap), domain)

similardomain(d::MappedDomain, ::Type{T}) where {T} =
    MappedDomain{T}(d.invmap, d.domain)

forward_map(d::MappedDomain) = rightinverse(d.invmap)
forward_map(d::MappedDomain, x) = rightinverse(d.invmap, x)

inverse_map(d::MappedDomain) = d.invmap
inverse_map(d::MappedDomain, y) = d.invmap(y)

"Map a domain with the inverse of the given map"
map_domain(map, domain::Domain) = _map_domain(map, domain)

# Fallback: we don't know anything about map, just try to invert
_map_domain(map, domain) = mapped_domain(inverse(map), domain)
_map_domain(map::Map{T}, domain::Domain{T}) where {T} =
    mapped_domain(inverse(map), domain)
# If map is a Map{T}, then verify and if necessary update T
function _map_domain(map::Map, domain)
    U = codomaintype(map, eltype(domain))
    if U == Union{}
        error("incompatible types of $(map) and $(domain)")
    end
    mapped_domain(inverse(convert(Map{U}, map)), convert(Domain{U}, domain))
end

==(a::MappedDomain, b::MappedDomain) = (a.invmap == b.invmap) && (superdomain(a) == superdomain(b))


"Make a mapped domain with the given inverse map"
mapped_domain(invmap, domain::Domain) = _mapped_domain(invmap, domain)

# We face the same task as in the constructor: attempt to identify T
# Here, we are more flexible, and attempt to do more conversions. We assume
# that users invoking the MappedDomain constructor know what they are doing,
# but users invoking mapped_domain just expect it to work.

# - we don't know anything about invmap, just pass it on
_mapped_domain(invmap, domain) = MappedDomain(invmap, domain)
# - invmap is a Map{T}: its codomaintype should match the eltype of the domain
# -- first, update the numtype
_mapped_domain(invmap::Map, domain) =
    _mapped_domain(invmap, domain, promote_type(numtype(invmap),numtype(domain)))
_mapped_domain(invmap::Map{T}, domain::Domain{S}, ::Type{U}) where {S,T,U} =
    _mapped_domain2(convert_numtype(invmap,U), convert_numtype(domain,U))
# -- then, ensure the codomaintype of the map equals the element type of the domain
_mapped_domain2(invmap, domain) = _mapped_domain2(invmap, domain, codomaintype(invmap), eltype(domain))
# --- it's okay
_mapped_domain2(invmap, domain, ::Type{T}, ::Type{T}) where {T} =
    MappedDomain(invmap, domain)
# --- it's not okay: attempt to convert the map (triggers e.g. when combining a scalar
#        LinearMap with a vector domain)
_mapped_domain2(invmap, domain, ::Type{S}, ::Type{T}) where {S,T} =
    MappedDomain(convert(Map{T}, invmap), domain)

# TODO: deprecate this syntax, it is confusing
(∘)(domain::Domain, invmap::Function) = mapped_domain(invmap, domain)
(∘)(domain::Domain, invmap::AbstractMap) = mapped_domain(invmap, domain)

# Avoid nested mapping domains, construct a composite map instead
# This assumes that the map types can be combined using \circ
mapped_domain(invmap, d::MappedDomain) = mapped_domain(inverse_map(d) ∘ invmap, superdomain(d))

boundary(d::MappedDomain) = _boundary(d, boundary(superdomain(d)), inverse_map(d))
_boundary(d::MappedDomain, superbnd, invmap) = MappedDomain(invmap, superbnd)
_boundary(d::MappedDomain, superbnd::UnionDomain, invmap) =
    UnionDomain(map(t->mapped_domain(invmap, t), components(superbnd)))

interior(d::MappedDomain) = _interior(d, superdomain(d), inverse_map(d))
_interior(d::MappedDomain, superdomain, invmap) = MappedDomain(invmap, interior(superdomain))
closure(d::MappedDomain) = _closure(d, superdomain(d), inverse_map(d))
_closure(d::MappedDomain, superdomain, invmap) = MappedDomain(invmap, closure(superdomain))

convert(::Type{MappedDomain}, d::Domain{T}) where {T} =
    MappedDomain{T}(mapto_canonical(d), canonicaldomain(d))
convert(::Type{MappedDomain{T}}, d::Domain) where {T} =
    MappedDomain{T}(mapto_canonical(d), canonicaldomain(d))


"A `ParametricDomain` stores the forward map of a mapped domain."
struct ParametricDomain{T,F,D} <: AbstractMappedDomain{T}
    fmap    ::  F
    domain  ::  D
end

ParametricDomain(fmap, domain::Domain) = ParametricDomain{codomaintype(fmap)}(fmap, domain)
ParametricDomain{T}(fmap, domain::Domain) where {T} =
    ParametricDomain{T,typeof(fmap),typeof(domain)}(fmap, domain)

similardomain(d::ParametricDomain, ::Type{T}) where {T} =
    ParametricDomain{T}(d.fmap, d.domain)

forward_map(d::ParametricDomain) = d.fmap
forward_map(d::ParametricDomain, x) = d.fmap(x)

function indomain(x, d::ParametricDomain)
    # To check for membership, we can't use the inverse map because it may not exist
    # We assume a left inverse exists, but the left inverse may be many-to-one.
    # So we also have to check whether the left-inverse-point maps back to x
    y = leftinverse(d.fmap, x)
    x2 = forward_map(d, y)
    isapprox(x, x2)
end

==(d1::ParametricDomain, d2::ParametricDomain) =
    forward_map(d1) == forward_map(d2) && superdomain(d1) == superdomain(d2)
hash(d::ParametricDomain, h::UInt) = hashrec(forward_map(d), superdomain(d), h)

"Return the domain that results from mapping the given domain."
parametric_domain(fmap, domain::Domain) = ParametricDomain(fmap, domain)
parametric_domain(fmap, domain::ParametricDomain) =
    parametric_domain(fmap ∘ forward_map(domain), superdomain(domain))

boundary(d::ParametricDomain) = _boundary(d, boundary(superdomain(d)), forward_map(d))
_boundary(d::ParametricDomain, superbnd, fmap) = ParametricDomain(fmap, superbnd)
_boundary(d::ParametricDomain, superbnd::UnionDomain, fmap) =
    UnionDomain(map(t -> parametric_domain(fmap, t), components(superbnd)))

interior(d::ParametricDomain) = _interior(d, superdomain(d), forward_map(d))
_interior(d::ParametricDomain, superdomain, fmap) = ParametricDomain(fmap, interior(superdomain))
closure(d::ParametricDomain) = _closure(d, superdomain(d), forward_map(d))
_closure(d::ParametricDomain, superdomain, fmap) = ParametricDomain(fmap, closure(superdomain))
