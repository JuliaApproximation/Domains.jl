
iscompact(d::TypedEndpointsInterval{:closed,:closed}) = true
iscompact(d::TypedEndpointsInterval) = false

isinterval(d::Domain) = false
isinterval(d::AbstractInterval) = true

hash(d::AbstractInterval, h::UInt) =
    hashrec(isleftopen(d), isrightopen(d), leftendpoint(d), rightendpoint(d), h)

Display.object_parentheses(::AbstractInterval) = true

approx_indomain(x, d::AbstractInterval, tolerance) =
    (x >= leftendpoint(d)-tolerance) && (x <= rightendpoint(d)+tolerance)

approx_indomain(x, d::TypedEndpointsInterval{:closed,:closed}, tolerance) =
    (x >= leftendpoint(d)-tolerance) && (x <= rightendpoint(d)+tolerance)
approx_indomain(x, d::TypedEndpointsInterval{:open,:closed}, tolerance) =
    (x > leftendpoint(d)-tolerance) && (x <= rightendpoint(d)+tolerance)
approx_indomain(x, d::TypedEndpointsInterval{:closed,:open}, tolerance) =
    (x >= leftendpoint(d)-tolerance) && (x < rightendpoint(d)+tolerance)
approx_indomain(x, d::TypedEndpointsInterval{:open,:open}, tolerance) =
    (x > leftendpoint(d)-tolerance) && (x < rightendpoint(d)+tolerance)


function point_in_domain(d::AbstractInterval)
    isempty(d) && throw(BoundsError())
    mean(d)
end

# For an interval of integers, try to find an integer point
function point_in_domain(d::AbstractInterval{T}) where {T<:Integer}
    isempty(d) && throw(BoundsError())
    x = round(T, mean(d))
    if x ∈ d
        x
    else
        # the mean is not inside the interval when rounded, this means
        # we have to choose an endpoint
        if isleftclosed(d)
            leftendpoint(d)
        elseif isrightclosed(d)
            rightendpoint(d)
        else
            # the interval is open and contains no integers in its interior
            throw(BoundsError())
        end
    end
end


isapprox(d1::AbstractInterval, d2::AbstractInterval) =
    isapprox(leftendpoint(d1), leftendpoint(d2); atol=default_tolerance(d1)) &&
    isapprox(rightendpoint(d1), rightendpoint(d2); atol=default_tolerance(d1))


boundary(d::AbstractInterval) = Point(leftendpoint(d)) ∪ Point(rightendpoint(d))
corners(d::AbstractInterval) = [leftendpoint(d), rightendpoint(d)]

normal(d::AbstractInterval, x) = (abs(minimum(d)-x) < abs(maximum(d)-x)) ? -one(eltype(d)) : one(eltype(d))

distance_to(d::AbstractInterval, x) = x ∈ d ? zero(eltype(d)) : min(abs(x-supremum(d)), abs(x-infimum(d)))

boundingbox(d::AbstractInterval) = d

volume(d::AbstractInterval) = width(d)

similar_interval(d::AbstractInterval, a, b) = similar_interval(d, promote(a, b)...)

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

closure(d::AbstractInterval) = ClosedInterval(endpoints(d)...)
closure(d::ClosedInterval) = d

interior(d::AbstractInterval) = OpenInterval(endpoints(d)...)

"""
Return an interval that is similar to the given interval, but with endpoints
`a` and `b` instead.
"""
similar_interval(d::ClosedFixedInterval{T}, a::S, b::S) where {S,T} =
    ClosedInterval{promote_type(float(T),S)}(a, b)

"The closed unit interval [0,1]."
struct UnitInterval{T} <: ClosedFixedInterval{T} end

UnitInterval() = UnitInterval{Float64}()

endpoints(d::UnitInterval{T}) where {T} = (zero(T), one(T))

similardomain(::UnitInterval, ::Type{T}) where {T} = UnitInterval{T}()


"The closed interval [-1,1]."
struct ChebyshevInterval{T} <: ClosedFixedInterval{T}
end

ChebyshevInterval() = ChebyshevInterval{Float64}()

endpoints(d::ChebyshevInterval{T}) where {T} = (-one(T),one(T))

similardomain(::ChebyshevInterval, ::Type{T}) where {T} = ChebyshevInterval{T}()

-(d::ChebyshevInterval) = d


"Map the interval [a,b] to the interval [c,d]."
interval_map(a, b, c, d) = AffineMap((d-c)/(b-a), c - a*(d-c)/(b-a))

canonicaldomain(d::ClosedInterval{T}) where {T} = ChebyshevInterval{float(T)}()
canonicaldomain(d::Interval{L,R,T}) where {L,R,T} = Interval{L,R,float(T)}(-1, 1)
canonicaldomain(d::FixedInterval) = d

mapfrom_canonical(d::Interval{L,R,T}) where {L,R,T} = interval_map(-one(T), one(T), endpoints(d)...)
mapfrom_canonical(d::FixedInterval{L,R,T}) where {L,R,T} = StaticIdentityMap{T}()

mapto(d1::D, d2::D) where {D <: FixedInterval} = identitymap(d1)
mapto(d1::D1, d2::D2) where {D1 <: FixedInterval,D2 <: FixedInterval} =
    interval_map(infimum(d1), supremum(d1), infimum(d2), supremum(d2))

"The half-open positive halfline `[0,∞)`."
struct HalfLine{T} <: FixedInterval{:closed,:open,T} end
HalfLine() = HalfLine{Float64}()

endpoints(d::HalfLine{T}) where {T} = (zero(T), T(Inf))

boundary(d::HalfLine) = Point(leftendpoint(d))

similardomain(::HalfLine, ::Type{T}) where {T} = HalfLine{T}()

# intercept and simplify the definition of IntervalSets
in(x, d::HalfLine) = x >= 0

approx_indomain(x, d::HalfLine, tolerance) = x >= -tolerance

function similar_interval(d::HalfLine{T}, a::S, b::S) where {T,S}
    @assert a == 0
    @assert isinf(b) && b > 0
    HalfLine{promote_type(float(T),S)}()
end

point_in_domain(d::HalfLine) = zero(eltype(d))


"The open negative halfline `(-∞,0)`."
struct NegativeHalfLine{T} <: FixedInterval{:open,:open,T} end
NegativeHalfLine() = NegativeHalfLine{Float64}()

similardomain(::NegativeHalfLine, ::Type{T}) where {T} = NegativeHalfLine{T}()

endpoints(d::NegativeHalfLine{T}) where {T} = (-T(Inf), zero(T))

boundary(d::NegativeHalfLine) = Point(rightendpoint(d))

# intercept and simplify the definition of IntervalSets
in(x, d::NegativeHalfLine) = x < 0

approx_indomain(x, d::NegativeHalfLine, tolerance) = x < tolerance

function similar_interval(d::NegativeHalfLine{T}, a::S, b::S) where {S,T}
    @assert isinf(a) && a < 0
    @assert b == 0
    NegativeHalfLine{promote_type(S,float(T))}()
end

point_in_domain(d::NegativeHalfLine) = -one(eltype(d))

similar_interval(d::Interval{L,R,T}, a::S, b::S) where {L,R,T,S} =
    Interval{L,R,promote_type(float(T),S)}(a, b)


#########################################
# A few set operations with known result
#########################################

# We define an exhaustive list of combinations of the four fixed intervals
# combined above in the routines 'intersect', 'union' and 'setdiff' where
# the output is known explicitly.

# Override the definition of Intervals.jl for FixedInterval's defined here
union(d1::FixedInterval, d2::FixedInterval) = uniondomain(d1, d2)
intersect(d1::FixedInterval, d2::FixedInterval) = intersectdomain(d1, d2)

uniondomain(d1::TypedEndpointsInterval, d2::TypedEndpointsInterval) =
    uniondomain(promote(d1,d2)...)

# type-unstable union of intervals. This function differs from `union` in
# IntervalSets.jl because that one throws an error if the intervals do not overlap
function uniondomain(d1::TypedEndpointsInterval{L1,R1,T}, d2::TypedEndpointsInterval{L2,R2,T}) where {L1,R1,L2,R2,T}
    a1 = leftendpoint(d1)
    b1 = rightendpoint(d1)
    a2 = leftendpoint(d2)
    b2 = rightendpoint(d2)

    # are they empty?
    isempty(d1) && return d2
    isempty(d2) && return d1
    # are they equal?
    d1 == d2 && return d1
    # does one lie within the other?
    a2 > a1 && b2 < b1 && return d1
    a1 > a2 && b1 < b2 && return d2
    # are they disjoint?
    b2 < a1 && return UnionDomain(d1, d2)   # return a UnionDomain for disjoint intervals
    b1 < a2 && return UnionDomain(d1, d2)
    # at this stage they must overlap
    if a1 < a2
        a = a1
        L = L1
    elseif a1 == a2
        a = a1
        L = (L1==L2==:open) ? :open : :closed
    else
        a = a2
        L = L2
    end
    if b1 > b2
        b = b1
        R = R1
    elseif b1 == b2
        b = b1
        R = (R1==R2==:open) ? :open : :closed
    else
        b = b2
        R = R2
    end
    Interval{L,R,T}(a, b)
end

function intersectdomain(d1::TypedEndpointsInterval, d2::TypedEndpointsInterval)
    # go back to the definition of IntervalSets.jl
    d = intersect(d1, d2)
    # but avoid returning an interval like 2..1
    isempty(d) ? EmptySpace{eltype(d)}() : d
end

# intersectdomain(d1::AbstractInterval, d2::AbstractInterval) = intersect(d1, d2)
# uniondomain(d1::AbstractInterval, d2::AbstractInterval) = union(d1, d2)

# Since fixed intervals are fully determined by their type,
# the result of intersect, union or setdiff is always known for two
# domains of the same type.
intersectdomain(d1::D, d2::D) where {D <: FixedInterval} = d1
uniondomain(d1::D, d2::D) where {D <: FixedInterval} = d1
setdiffdomain(d1::D, d2::D) where {D <: FixedInterval} = EmptySpace{eltype(D)}()

# [0,1] ∩ [-1,1] = [0,1]
intersectdomain(d1::UnitInterval{T}, d2::ChebyshevInterval{T}) where {T} = UnitInterval{T}()
intersectdomain(d1::ChebyshevInterval{T}, d2::UnitInterval{T}) where {T} = UnitInterval{T}()
# [0,1] ∩ [0,∞) = [0,1]
intersectdomain(d1::UnitInterval{T}, d2::HalfLine{T}) where {T} = UnitInterval{T}()
intersectdomain(d1::HalfLine{T}, d2::UnitInterval{T}) where {T} = UnitInterval{T}()
# [0,1] ∩ (-∞,0) = {}
intersectdomain(d1::UnitInterval{T}, d2::NegativeHalfLine{T}) where {T} = EmptySpace{T}()
intersectdomain(d1::NegativeHalfLine{T}, d2::UnitInterval{T}) where {T} = EmptySpace{T}()
# [-1,1] ∩ [0,∞) = [0,1]
intersectdomain(d1::ChebyshevInterval{T}, d2::HalfLine{T}) where {T} = UnitInterval{T}()
intersectdomain(d1::HalfLine{T}, d2::ChebyshevInterval{T}) where {T} = UnitInterval{T}()
# [0,∞) ∩ (-∞,0) = {}
intersectdomain(d1::HalfLine{T}, d2::NegativeHalfLine{T}) where {T} = EmptySpace{T}()
intersectdomain(d1::NegativeHalfLine{T}, d2::HalfLine{T}) where {T} = EmptySpace{T}()

# [0,1] ∪ [-1,1] = [-1,1]
uniondomain(d1::UnitInterval{T}, d2::ChebyshevInterval{T}) where {T} = ChebyshevInterval{T}()
uniondomain(d1::ChebyshevInterval{T}, d2::UnitInterval{T}) where {T} = ChebyshevInterval{T}()
# [0,1] ∪ [0,∞) = [0,∞)
uniondomain(d1::UnitInterval{T}, d2::HalfLine{T}) where {T} = HalfLine{T}()
uniondomain(d1::HalfLine{T}, d2::UnitInterval{T}) where {T} = HalfLine{T}()

# (-∞,0) ∪ [0,∞) = (-∞,∞)
# Note: T<:real to ensure that FullSpace{T} is not larger than intended.
uniondomain(d1::NegativeHalfLine{T}, d2::HalfLine{T}) where {T<:Real} = FullSpace{T}()
uniondomain(d1::HalfLine{T}, d2::NegativeHalfLine{T}) where {T<:Real} = FullSpace{T}()


# [0,1] ∖ [-1,1] = {}
setdiffdomain(d1::UnitInterval{T}, d2::ChebyshevInterval{T}) where {T} = EmptySpace{T}()
# [0,1] ∖ [0,∞) = {}
setdiffdomain(d1::UnitInterval{T}, d2::HalfLine{T}) where {T} = EmptySpace{T}()
# [0,1] ∖ (-∞,0) = [0,1]
setdiffdomain(d1::UnitInterval{T}, d2::NegativeHalfLine{T}) where {T} = UnitInterval{T}()
# [-1,1] ∖ (-∞,0) = [0,1]
setdiffdomain(d1::ChebyshevInterval{T}, d2::NegativeHalfLine{T}) where {T} = UnitInterval{T}()
# [0,∞) ∖ (-∞,0) = [0,∞)
setdiffdomain(d1::HalfLine{T}, d2::NegativeHalfLine{T}) where {T} = HalfLine{T}()
# (-∞,0) ∖ [0,1] = (-∞,0)
setdiffdomain(d1::NegativeHalfLine{T}, d2::UnitInterval{T}) where {T} = NegativeHalfLine{T}()
# (-∞,0) ∖ [0,∞) = (-∞,0)
setdiffdomain(d1::NegativeHalfLine{T}, d2::HalfLine{T}) where {T} = NegativeHalfLine{T}()


#################################
# Conversions between intervals
#################################
# Avoid depcrecated warning: Warning: Constructors no longer fall back to `convert`.
# example: A constructor `AbstractInterval{Float64}(::IntervalSets.ClosedInterval{Float64})` should be defined instead.

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


function show(io::IO, d::ChebyshevInterval)
    print(io, Interval(d))
    print(io, " (Chebyshev)")
end
function show(io::IO, d::UnitInterval)
    print(io, Interval(d))
    print(io, " (Unit)")
end
function show(io::IO, d::HalfLine)
    print(io, Interval(d))
    print(io, " (HalfLine)")
end


########################
# Arithmetic operations
########################

function setdiffdomain(d1::AbstractInterval{T}, d2::AbstractInterval{T}) where T
    a1 = leftendpoint(d1)
    b1 = rightendpoint(d1)
    a2 = leftendpoint(d2)
    b2 = rightendpoint(d2)

    isempty(d1) && return d1
    isempty(d2) && return d1
    b1 < a2 && return d1
    a1 < a2 ≤ b1 ≤ b2 && return (a1..a2)
    a1 < a2 ≤ b2 < b1 && return uniondomain(a1..a2, b2..b1)
    a2 ≤ a1 < b2 < b1 && return (b2..b1)
    a2 ≤ a1 ≤ b1 ≤ b2 && return EmptySpace{T}()

    @assert b2 ≤ a1
    d1
end

setdiffdomain(d1::AbstractInterval, d2::AbstractInterval) =
    setdiffdomain(promote(d1,d2)...)

switch_open_closed(d::AbstractInterval) = d
switch_open_closed(d::Interval{L,R,T}) where {L,R,T} =
    Interval{R,L,T}(leftendpoint(d),rightendpoint(d))

function map_domain(map::AbstractAffineMap{<:Number}, domain::AbstractInterval)
    le = map(leftendpoint(domain))
    re = map(rightendpoint(domain))
    if le<re
        similar_interval(domain,le,re)
    else
        similar_interval(switch_open_closed(domain),re,le)
    end
end

mapped_domain(invmap::AbstractAffineMap{<:Number}, domain::AbstractInterval) =
    map_domain(inverse(invmap), domain)

# Preserve maps when acting on the half line and negative halfline, because
# the map may not always preserve the domain type
map_domain(map::AbstractAffineMap{<:Number}, domain::HalfLine) = MappedDomain(inverse(map), domain)
map_domain(map::AbstractAffineMap{<:Number}, domain::NegativeHalfLine) = MappedDomain(inverse(map), domain)
