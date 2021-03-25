##################
### Segment
##################


"""
	Segment(a,b)

represents a line segment from `a` to `b`.  In the case where `a` and `b`
are real and `a < b`, then this is is equivalent to an `Interval(a,b)`.
"""
struct Segment{T} <: Domain{T}
	a::T
	b::T
	Segment{T}(a,b) where {T} = new{T}(a,b)
end

const IntervalOrSegment{T} = Union{AbstractInterval{T}, Segment{T}}

Segment(a::Complex{IT1}, b::Complex{IT2}) where {IT1<:Integer,IT2<:Integer} =
	Segment(ComplexF64(a), ComplexF64(b)) #convenience method
Segment(a::Integer, b::Integer) = Segment(Float64(a),Float64(b)) #convenience method
Segment(a::Complex{IT}, b) where {IT<:Integer} = Segment(ComplexF64(a),b) #convenience method
Segment(a, b::Complex{IT}) where {IT<:Integer} = Segment(a,ComplexF64(b)) #convenience method
Segment(a, b) = Segment{promote_type(typeof(a),typeof(b))}(a,b)
Segment(a::Tuple, b::Tuple) = Segment(Vec(a...),Vec(b...))

convert(::Type{Domain{T}}, d::Segment) where {T<:Number} = Segment{T}(leftendpoint(d),rightendpoint(d))
convert(::Type{Domain{T}}, d::Segment) where {T<:SVector} = Segment{T}(leftendpoint(d),rightendpoint(d))
convert(::Type{Segment{T}}, d::Segment) where {T<:Number} = Segment{T}(leftendpoint(d),rightendpoint(d))
convert(::Type{Segment}, d::AbstractInterval) = Segment(leftendpoint(d), rightendpoint(d))
convert(::Type{Segment{T}}, d::AbstractInterval) where T =convert(Segment{T}, convert(Segment, d))
convert(::Type{Interval}, d::Segment{<:Real}) = d.a < d.b ? d.a .. d.b : d.b .. d.a

Segment(d::AbstractInterval) = convert(Segment, d)
Interval(d::Segment) = convert(Interval, d)

@inline leftendpoint(d::Segment) = d.a
@inline rightendpoint(d::Segment) = d.b
@inline endpoints(d::Segment) = d.a, d.b

==(d::Segment, m::Segment) = leftendpoint(d) == leftendpoint(m) && rightendpoint(d) == rightendpoint(m)
function isapprox(d::Segment, m::Segment)
    tol=10E-12
    norm(leftendpoint(d)-leftendpoint(m))<tol && norm(rightendpoint(d)-rightendpoint(m))<tol
end

for op in (:(==), :isapprox)
    @eval begin
        $op(d::Segment, m::AbstractInterval) = $op(d, Segment(m))
        $op(m::AbstractInterval, d::Segment) = $op(Segment(m), d)
    end
end

@inline minimum(d::Segment) = min(leftendpoint(d),rightendpoint(d))
@inline maximum(d::Segment) = max(leftendpoint(d),rightendpoint(d))

isempty(d::Segment) = isapprox(leftendpoint(d), rightendpoint(d); atol=200eps(eltype(d)))

issubset(a::Segment,b::Segment) = leftendpoint(a) ∈ b && rightendpoint(a) ∈ b

arclength(d::AbstractInterval) = width(d)
arclength(d::Segment) = norm(complexlength(d))
complexlength(d::IntervalOrSegment) = rightendpoint(d)-leftendpoint(d)
mean(d::IntervalOrSegment) = (rightendpoint(d)+leftendpoint(d))/2
angle(d::IntervalOrSegment) = angle(complexlength(d))
sign(d::IntervalOrSegment) = sign(complexlength(d))

indomain(x, d::Segment{<:Real}) = indomain(x, Interval(d))
function indomain(x, d::Segment)
    xda = Segment(leftendpoint(d), x)
    angle(xda) == angle(d) && arclength(xda) ≤ arclength(d)
end

intersect(a::Segment{<:Real}, b::Segment{<:Real}) = intersect(Interval(a), Interval(b))
intersect(a::AbstractInterval, b::Segment{<:Real}) = intersect(a, Interval(b))
intersect(a::Segment{<:Real}, b::AbstractInterval) = intersect(Interval(a), b)
setdiff(a::Segment{<:Real}, b::Segment{<:Real})  = setdiff(Interval(a), Interval(b))
setdiff(a::AbstractInterval, b::Segment{<:Real})  = setdiff(a, Interval(b))
setdiff(a::Segment{<:Real}, b::AbstractInterval)  = setdiff(Interval(a), b)

isless(d1::Segment{T1},d2::Segment{T2}) where {T1<:Real,T2<:Real} =
    d1 ≤ leftendpoint(d2) && d1 ≤ rightendpoint(d2)
isless(d1::Segment{T},x::Real) where {T<:Real} = leftendpoint(d1) ≤ x && rightendpoint(d1) ≤ x
isless(x::Real,d1::Segment{T}) where {T<:Real} = x ≤ leftendpoint(d1) && x ≤ rightendpoint(d1)
