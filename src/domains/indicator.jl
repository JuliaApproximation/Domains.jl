
"""
Supertype of domains that are defined by an indicator function.

An indicator function is a function `f : S -> [0,1]` that indicates membership
of `x` to a domain `D` with `D ⊂ S`. The indicator function corresponds exactly
to the `in` function of a domain: `f(x) = x ∈ D`.

Concrete subtypes of `AbstractIndicatorFunction` store a representation of this
indicator function and implement `in` using that representation, rather than
implementing `in` directly.
"""
abstract type AbstractIndicatorFunction{T} <: Domain{T} end

"The indicator function of a domain is the function `f(x) = x ∈ D`."
indicatorfunction(d::Domain) = x -> x ∈ d

indomain(x, d::AbstractIndicatorFunction) = _indomain(x, d, indicatorfunction(d))
_indomain(x, d::AbstractIndicatorFunction, f) = f(x)

show(io::IO, d::AbstractIndicatorFunction) =
    print(io, "indicator domain defined by function f = $(indicatorfunction(d))")


"An `IndicatorFunction` is a domain that implements `f(x) = x ∈ D` by storing `f`."
struct IndicatorFunction{T,F} <: AbstractIndicatorFunction{T}
    f       ::  F
end

IndicatorFunction(f) = IndicatorFunction{Float64}(f)
IndicatorFunction{T}(f::F) where {T,F} = IndicatorFunction{T,F}(f)

indicatorfunction(d::IndicatorFunction) = d.f

similardomain(d::IndicatorFunction, ::Type{T}) where {T} = IndicatorFunction{T}(d.f)

convert(::Type{IndicatorFunction}, d::AbstractIndicatorFunction) = d
convert(::Type{IndicatorFunction}, d::Domain{T}) where {T} =
    IndicatorFunction{T}(indicatorfunction(d))

==(d1::IndicatorFunction, d2::IndicatorFunction) = indicatorfunction(d1)==indicatorfunction(d2)

intersectdomain1(d1::IndicatorFunction, d2) = BoundedIndicatorFunction(d1.f, d2)
intersectdomain2(d1, d2::IndicatorFunction) = BoundedIndicatorFunction(d2.f, d1)

"An indicator function with a known bounding domain."
struct BoundedIndicatorFunction{F,D,T} <: AbstractIndicatorFunction{T}
    f       ::  F
    domain  ::  D
end

BoundedIndicatorFunction(f::F, domain::D) where {F,T,D<:Domain{T}} =
    BoundedIndicatorFunction{F,D,T}(f, domain)

indicatorfunction(d::BoundedIndicatorFunction) = d.f

boundingdomain(d::BoundedIndicatorFunction) = d.domain

indomain(x, d::BoundedIndicatorFunction) = in(x, boundingdomain(d)) && d.f(x)

==(d1::BoundedIndicatorFunction, d2::BoundedIndicatorFunction) =
    indicatorfunction(d1)==indicatorfunction(d2) && boundingdomain(d1)==boundingdomain(d2)
hash(d::BoundedIndicatorFunction, h::UInt) =
    hashrec(indicatorfunction(d), boundingdomain(d), h)

similardomain(d::BoundedIndicatorFunction, ::Type{T}) where {T} =
    BoundedIndicatorFunction(d.f, convert(Domain{T}, d.domain))

Domain(gen::Base.Generator) = generator_domain(gen)

generator_domain(gen::Base.Generator{<:Domain}) = BoundedIndicatorFunction(gen.f, gen.iter)
generator_domain(gen::Base.Generator{<:Base.Iterators.ProductIterator}) =
    productgenerator_domain(gen, gen.iter.iterators)

function productgenerator_domain(gen, domains::Tuple{Vararg{Domain,N} where N})
    domain = TupleProductDomain(gen.iter.iterators)
    BoundedIndicatorFunction(gen.f, domain)
end

boundingbox(d::BoundedIndicatorFunction) = boundingbox(boundingdomain(d))

function show(io::IO, d::BoundedIndicatorFunction)
    print(io, "indicator function bounded by: ")
    show(io, boundingdomain(d))
end
