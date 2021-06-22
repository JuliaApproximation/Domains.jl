
"""
    canonicaldomain([ctype::CanonicalType, ]d::Domain)

Return an associated canonical domain, if any, of the given domain.

Optionally, a canonical type argument may specify an alternative canonical domain.
Canonical domains help with establishing equality between domains, with finding
maps between domains and with finding parameterizations.

For example, the canonical domain of an Interval `[a,b]` is the interval `[-1,1]`.
"""
canonicaldomain(d::Domain) = d

"""
Check whether the given arguments are different. They are different when
they have a different type, or if they have the same type but they are not equal.
"""
isdifferentfrom(d1::D, d2::D) where D = !(d1==d2)
isdifferentfrom(d1, d2) = true

"Does the domain have a canonical domain?"
hascanonicaldomain(d) = isdifferentfrom(d, canonicaldomain(d))

identitymap(d) = IdentityMap{eltype(d)}(dimension(d))

"Return a map to a domain from its canonical domain."
mapfrom_canonical(d) = identitymap(d)
mapfrom_canonical(d, x) = mapfrom_canonical(d)(x)

"Return a map from the domain to its canonical domain."
mapto_canonical(d) = leftinverse(mapfrom_canonical(d))
mapto_canonical(d, x) = mapto_canonical(d)(x)


"Supertype of kinds of canonical domains."
abstract type CanonicalType end

canonicaldomain(ctype::CanonicalType, d) = d
hascanonicaldomain(ctype::CanonicalType, d) = isdifferentfrom(d, canonicaldomain(ctype, d))

mapfrom_canonical(ctype::CanonicalType, d, x) = mapfrom_canonical(ctype, d)(x)
mapto_canonical(ctype::CanonicalType, d, x) = mapto_canonical(ctype, d)(x)


"A canonical domain that is equal but simpler (e.g. a 1-dimensional ball is an interval)."
struct Equal <: CanonicalType end

canonicaldomain(::Equal, d) = d
mapfrom_canonical(::Equal, d) = identitymap(d)
mapto_canonical(::Equal, d) = leftinverse(mapto_canonical(Equal(), d))

simplify(d) = canonicaldomain(Equal(), d)
simplifies(d) = hascanonicaldomain(Equal(), d)

"A canonical domain that is isomorphic but may have different element type."
struct Isomorphic <: CanonicalType end

canonicaldomain(::Isomorphic, d) = canonicaldomain(Equal(), d)
mapfrom_canonical(::Isomorphic, d) = mapfrom_canonical(Equal(), d)
mapto_canonical(::Isomorphic, d) = leftinverse(mapto_canonical(Isomorphic(), d))

canonicaldomain(::Isomorphic, d::Domain{SVector{1,T}}) where {T} =
    convert(Domain{T}, d)
mapfrom_canonical(::Isomorphic, d::Domain{SVector{1,T}}) where {T} =
    NumberToVector{T}()

canonicaldomain(::Isomorphic, d::Domain{NTuple{N,T}}) where {N,T} =
    convert(Domain{SVector{N,T}}, d)
mapfrom_canonical(::Isomorphic, d::Domain{NTuple{N,T}}) where {N,T} =
    VectorToTuple{N,T}()


"A parameter domain that can be mapped to the domain."
struct Parameterization <: CanonicalType end

canonicaldomain(ctype::Parameterization, d) =
    hascanonicaldomain(d) ? canonicaldomain(ctype, canonicaldomain(d)) : d
mapfrom_canonical(ctype::Parameterization, d) =
    hascanonicaldomain(d) ? mapfrom_canonical(d) ∘ mapfrom_canonical(ctype, canonicaldomain(d)) : mapfrom_canonical(d)
mapto_canonical(ctype::Parameterization, d) = leftinverse(mapfrom_canonical(ctype, d))

# We define some convenience functions:
"Return a parameter domain which supports a `parameterization`."
parameterdomain(d::Domain) = canonicaldomain(Parameterization(), d)

"Return a parameterization of the given domain."
parameterization(d::Domain) = mapfrom_canonical(Parameterization(), d)

"Does the domain have a parameterization?"
hasparameterization(d) = hascanonicaldomain(Parameterization(), d)

mapfrom_parameterdomain(d::Domain) = mapfrom_canonical(Parameterization(), d)
mapto_parameterdomain(d::Domain) = mapto_canonical(Parameterization(), d)


"Return a map from domain `d1` to domain `d2`."
mapto(d1, d2) = mapto1(d1, d2)
mapto(d1::D, d2::D) where {D} = d1 == d2 ? identitymap(d1) : mapto1(d1,d2)

# simplify the first argument
mapto1(d1, d2) =
    hasparameterization(d1) ? mapto(parameterdomain(d1), d2) ∘ mapto_parameterdomain(d1) : mapto2(d1, d2)
# simplify the second argument
mapto2(d1, d2) =
    hasparameterization(d2) ? mapfrom_parameterdomain(d2) ∘ mapto(d1, parameterdomain(d2)) : no_known_mapto(d1,d2)

no_known_mapto(d1, d2) = d1 == d2 ? identitymap(d1) : error("No map known between $(d1) and $(d2).")



## Equality for domains

==(d1::Domain, d2::Domain) = isequal1(d1, d2)
# simplify the first argument
isequal1(d1, d2) = simplifies(d1) ? simplify(d1)==d2 : isequal2(d1, d2)
# simplify the second argument
isequal2(d1, d2) = simplifies(d2) ? d1==simplify(d2) : false
