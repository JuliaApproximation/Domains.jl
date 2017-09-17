# simple.jl
# A collection of simple domains.

###########################
# The unit ball and sphere
###########################

"""
The unit ball (of radius 1) in `N` dimensions.
"""
struct UnitBall{N,T} <: EuclideanDomain{N,T}
end

const Disk{T} = UnitBall{2,T}
const Ball{T} = UnitBall{3,T}

indomain(x, ::UnitBall) = norm(x) <= 1

disk(::Type{T} = Float64) where {T} = Disk{T}()
disk(radius::Number) = radius * disk(typeof(radius))
disk(radius::Number, center::AbstractVector) = disk(radius) + center

ball(::Type{T} = Float64) where {T} = Ball{T}()
ball(radius::Number) = radius * ball(typeof(radius))
ball(radius::Number, center::AbstractVector) = ball(radius) + center

show(io::IO, d::UnitBall{N}) where {N} = print(io, "the $(N)-dimensional unit ball")




###########################
# An n-dimensional simplex
###########################

struct UnitSimplex{N,T} <: EuclideanDomain{N,T}
end

indomain(x, ::UnitSimplex) = mapreduce( t-> t >= 0, &, x) && norm(x,1) <= 1

simplex(::Type{Val{N}}, ::Type{T} = Float64) where {T,N} = UnitSimplex{N,T}()


#########################
# An n-dimensional cube
#########################

cube(::Type{Val{N}}, ::Type{T} = Float64) where {N,T} = cartesianproduct(UnitInterval{T}(), Val{N})

cube() = cube(Val{3})

rectangle(a, b, c, d) = interval(a,b) × interval(c,d)

cube(a, b, c, d, e, f) = interval(a,b) × interval(c,d) × interval(e,f)

# This one is not type-stable
cube(a::NTuple{N,T}, b::NTuple{N,T}) where {N,T} = ProductDomain(map((ai,bi)->interval(T,ai,bi), a, b)...)
# This one isn't either
cube(a::AbstractVector{T}, b::AbstractVector{T}) where {T} = cube(tuple(a...), tuple(b...))

# const Square{T} = UnitCube{2,T}
# const Cube{T} = UnitCube{3,T}}



##############
# A cylinder
##############

cylinder(::Type{T} = Float64) where {T} = disk(T) × unitinterval(T)

cylinder(radius::T, length::T) where {T} = disk(radius) × interval(0,length)

doc"""
    Point(x)

represents a single point at `x`.
"""
struct Point{T} <: Domain{T}
    x::T
end

convert(::Type{Number}, d::Point) = d.x
convert(::Type{<:Number}, d::Point) = N(d.x)

convert(::Type{Domain}, c::Number) = Point(c)
convert(::Type{Domain{T}}, c::Number) where T = Point(convert(T,c))



==(a::Point,b::Point) = a.x == b.x
indomain(x, d::Point) = x == d.x

for op in (:*,:+,:-)
    @eval begin
        $op(c::Number, d::Point)  = Point($op(c,d.x))
        $op(d::Point,  c::Number) = Point($op(d.x,c))
    end
end


/(d::Point,  c::Number) = Point(d.x/c)
\(c::Number, d::Point)  = Point(c\d.x)

for op in (:+,:-)
    @eval $op(a::Point, b::Point) = Point($op(a.x,b.x))
end


for op in (:*,:+)
    @eval begin
        $op(a::Point,v::AbstractVector) = map(y->$op(a,y),v)
        $op(v::AbstractVector,a::Point) = map(y->$op(y,a),v)
    end
end
