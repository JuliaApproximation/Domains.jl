# interval.jl


"A half-open interval `(a,b]`."
const HalfOpenLeftInterval{T} = Interval{:open,:closed,T}

"A half-open interval `[a,b)`."
const HalfOpenRightInterval{T} = Interval{:closed,:open,T}


iscompact(d::TypedEndpointsInterval{:closed,:closed}) = true
iscompact(d::TypedEndpointsInterval) = false

indomain(x, d::AbstractInterval) = x ∈ d


##################
### An interval
##################

approx_indomain(x, d::AbstractInterval, tolerance) =
    (x <= rightendpoint(d)+tolerance) && (x >= leftendpoint(d)-tolerance)


function point_in_domain(d::AbstractInterval)
    isempty(d) && throw(BoundsError())
    mean(d)
end

isapprox(d1::AbstractInterval, d2::AbstractInterval; kwds...) =
    isapprox(leftendpoint(d1), leftendpoint(d2); kwds...) &&
    isapprox(rightendpoint(d1), rightendpoint(d2); kwds...)


boundary(d::TypedEndpointsInterval{:closed,:closed}) = Point(leftendpoint(d)) ∪ Point(rightendpoint(d))
boundary(d::TypedEndpointsInterval{:closed,:open}) = Point(leftendpoint(d))
boundary(d::TypedEndpointsInterval{:open,:closed}) = Point(rightendpoint(d))
boundary(d::TypedEndpointsInterval{:open,:open,T}) where T = EmptySpace{T}()
boundary(d::AbstractInterval) = boundary(Interval(d))


# We extend some functionality of intervals to mapped intervals
const MappedInterval{D <: AbstractInterval,T} = MappedDomain{D,T}

endpoints(d::MappedInterval) = forward_map(d) .* endpoints(source(d))
for op in (:leftendpoint, :rightendpoint)
    @eval $op(d::MappedInterval) = forward_map(d) * $op(source(d))
end

for op in (:isleftopen, :isrightopen, :isleftclosed, :isrightclosed)
    @eval $op(d::MappedInterval) = $op(source(d))
end


## Some special intervals follow, e.g.:
# - the unit interval [0,1]
# - the 'Chebyshev' interval [-1,1]
# - ...
# Unlike a generic interval, these specific intervals have no data, and only one
# type parameter (T).

"""
The abstract type `FixedInterval` is the supertype of intervals with endpoints
determined by the type, rather than field values. Examples include `UnitInterval`
and `ChebyshevInterval`.
"""
abstract type FixedInterval{L,R,T} <: TypedEndpointsInterval{L,R,T} end
const ClosedFixedInterval{T} = FixedInterval{:closed,:closed,T}

"""
Return an interval that is similar to the given interval, but with endpoints
`a` and `b` instead.
"""# Assume a closed interval by default
similar_interval(d::ClosedFixedInterval{T}, a, b) where {T} = ClosedInterval{float(T)}(a, b)


"The closed unit interval [0,1]."
struct UnitInterval{T} <: ClosedFixedInterval{T} end

UnitInterval() = UnitInterval{Float64}()
unitinterval(::Type{T} = Float64) where {T} = UnitInterval{T}()

endpoints(d::UnitInterval{T}) where {T} = (zero(T), one(T))


"The closed interval [-1,1]."
struct ChebyshevInterval{T} <: ClosedFixedInterval{T}
end

ChebyshevInterval() = ChebyshevInterval{Float64}()

endpoints(d::ChebyshevInterval{T}) where {T} = (-one(T),one(T))

real_line(::Type{T} = Float64) where {T <: AbstractFloat} = FullSpace{T}()


"The half-open positive halfline `[0,∞)`."
struct Halfline{T} <: FixedInterval{:closed,:open,T}
end

halfline(::Type{T} = Float64) where {T <: AbstractFloat} = Halfline{T}()

endpoints(d::Halfline{T}) where {T} = (zero(T), T(Inf))


indomain(x, d::Halfline) = x >= 0

function similar_interval(d::Halfline, a, b)
    @assert a == 0
    @assert isinf(b) && b > 0
    d
end

point_in_domain(d::Halfline) = zero(eltype(d))


"The open negative halfline `(-∞,0)`."
struct NegativeHalfline{T} <: FixedInterval{:open,:open,T}
end

negative_halfline(::Type{T} = Float64) where {T <: AbstractFloat} = NegativeHalfline{T}()

endpoints(d::NegativeHalfline{T}) where {T} = (-T(Inf), zero(T))


# Open at both endpoints


indomain(x, d::NegativeHalfline) = x < 0

function similar_interval(d::NegativeHalfline, a, b)
    @assert isinf(a) && a < 0
    @assert b == 0
    d
end

point_in_domain(d::NegativeHalfline) = -one(eltype(d))


similar_interval(d::Interval{L,R,T}, a, b) where {L,R,T} =
    Interval{L,R,float(T)}(a, b)


#################################
# Conversions between intervals
#################################
# Avoid depcrecated warning: Warning: Constructors no longer fall back to `convert`.
# example: A constructor `Domains.AbstractInterval{Float64}(::IntervalSets.ClosedInterval{Float64})` should be defined instead.

function convert(::Type{UnitInterval{T}}, d::AbstractInterval) where {T}
    endpoints(d) == (0,1) || throw(InexactError(:convert,UnitInterval,d))
    UnitInterval{T}()
end

function convert(::Type{ChebyshevInterval{T}}, d::AbstractInterval) where {T}
    endpoints(d) == (-1,1)|| throw(InexactError(:convert,ChebyshevInterval,d))
    ChebyshevInterval{T}()
end

UnitInterval{T}(d::AbstractInterval) where T = convert(UnitInterval{T}, d)
ChebyshevInterval{T}(d::AbstractInterval) where T = convert(ChebyshevInterval{T}, d)


########################
# Arithmetic operations
########################

# Some computations with intervals simplify without having to use a mapped domain.
# This is only the case for Interval{L,R,T}, and not for any of the FixedIntervals
# because the endpoints of the latter are, well, fixed.

-(d::ChebyshevInterval) = d
-(d::AbstractInterval) = similar_interval(d, -rightendpoint(d), -leftendpoint(d))

for op in (:+, :-)
    @eval $op(d::AbstractInterval, x::Real) = similar_interval(d, $op(leftendpoint(d),x), $op(rightendpoint(d),x))
end

+(x::Real, d::AbstractInterval) = similar_interval(d, x+leftendpoint(d), x+rightendpoint(d))
-(x::Real, d::AbstractInterval) = similar_interval(d, x-rightendpoint(d), x-leftendpoint(d))

for op in (:*, :/)
    @eval function $op(d::AbstractInterval, x::Real)
        if x ≥ 0 # -{x : 0 ≤ x ≤ 1} should be {x : -1 ≤ x ≤ 0}, not empty set {x : 0 ≤ x ≤ -1}
            similar_interval(d, $op(leftendpoint(d),x), $op(rightendpoint(d),x))
        else
            similar_interval(d, $op(rightendpoint(d),x), $op(leftendpoint(d),x))
        end
    end
end

for op in (:*, :\)
    @eval function $op(x::Real, d::AbstractInterval)
        if x ≥ 0 # -{x : 0 ≤ x ≤ 1} should be {x : -1 ≤ x ≤ 0}, not empty set {x : 0 ≤ x ≤ -1}
            similar_interval(d, $op(x,leftendpoint(d)), $op(x,rightendpoint(d)))
        else
            similar_interval(d, $op(x,rightendpoint(d)), $op(x,leftendpoint(d)))
        end
    end
end


show(io::IO, d::AbstractInterval) = print(io, "an interval of type ", typeof(d), " ",
                                            isleftopen(d) ? "(" : "[", leftendpoint(d), ", ",
                                            rightendpoint(d), isrightopen(d) ? ")" : "]")
show(io::IO, d::ChebyshevInterval) = print(io,  "the Chebyshev interval [", leftendpoint(d), ", ", rightendpoint(d), "]")
show(io::IO, d::UnitInterval) = print(io,  "the unit interval [", leftendpoint(d), ", ", rightendpoint(d), "]")
show(io::IO, d::Halfline) = print(io,  "the half-line [", leftendpoint(d), ", ", rightendpoint(d), ")")

function setdiff(d1::AbstractInterval{T}, d2::AbstractInterval{T}) where T
    a1 = leftendpoint(d1)
    b1 = rightendpoint(d1)
    a2 = leftendpoint(d2)
    b2 = rightendpoint(d2)

    isempty(d1) && return d1
    isempty(d2) && return d1
    b1 < a2 && return d1
    a1 < a2 ≤ b1 ≤ b2 && return (a1 .. a2)
    a1 < a2 ≤ b2 < b1 && return UnionDomain(a1 .. a2) ∪ UnionDomain(b2 .. b1)
    a2 ≤ a1 < b2 < b1 && return (b2 .. b1)
    a2 ≤ a1 ≤ b1 ≤ b2 && return EmptySpace{T}()

    @assert b2 ≤ a1
    d1
end

setdiff(d1::AbstractInterval, d2::AbstractInterval) = setdiff(promote(d1,d2)...)

function *(map::AffineMap, domain::AbstractInterval)
    le = map*leftendpoint(domain)
    re = map*rightendpoint(domain)
    if le<re
        similar_interval(domain,le,re)
    else
        similar_interval(domain,re,le)
    end
end
