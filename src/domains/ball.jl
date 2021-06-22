
# The type hierarchy is as follows:
# abstract Ball
# |-> abstract UnitBall: radius is 1
#     |-> StaticUnitBall: dimension is part of type
#     |-> DynamicUnitBall: dimension is specified by int field
# |-> GenericBall: stores center and radius. Here, the dimension can be
#     obtained from the center, so there is only one type.
# There are aliases for SVector{N,T} and Vector{T}:
#   EuclideanBall, VectorBall,
#   EuclideanUnitBall (StaticUnitBall), VectorUnitBall (DynamicUnitBall).

"""
Supertype of balls for which elements satisfy `norm(x) < radius(ball)` (open ball)
or `norm(x) <= radius(ball)` (closed ball).
"""
abstract type Ball{T,C} <: Domain{T} end

"A ball in a fixed N-dimensional Euclidean space."
const EuclideanBall{N,T,C} = Ball{SVector{N,T},C}

"A ball with vector elements of variable length."
const VectorBall{T,C} = Ball{Vector{T},C}

iscompatiblepair(x::AbstractVector, d::VectorBall) = length(x) == dimension(d)

const ClosedBall{T} = Ball{T,:closed}
const OpenBall{T} = Ball{T,:open}

isclosedset(::ClosedBall) = true
isclosedset(::OpenBall) = false
isopenset(d::Ball) = !isclosedset(d)

# Convenience constructors for the abstract type
Ball() = UnitBall()
Ball{T}() where {T} = UnitBall{T}()
Ball{T,C}() where {T,C} = UnitBall{T,C}()
Ball(radius::Number) = GenericBall(radius)
Ball{T}(radius::Number) where {T} = GenericBall{T}(radius)
Ball{T,C}(radius::Number) where {T,C} = GenericBall{T,C}(radius)
Ball(radius::Number, center) = GenericBall(radius, center)
Ball{T}(radius::Number, center) where {T} = GenericBall{T}(radius, center)
Ball{T,C}(radius::Number, center) where {T,C} = GenericBall{T,C}(radius, center)


indomain(x, d::OpenBall) = norm(x-center(d)) < radius(d)
indomain(x, d::ClosedBall) = norm(x-center(d)) <= radius(d)
approx_indomain(x, d::OpenBall, tolerance) = norm(x-center(d)) < radius(d)+tolerance
approx_indomain(x, d::ClosedBall, tolerance) = norm(x-center(d)) <= radius(d)+tolerance

# A closed ball always contains at least its center point.
isempty(d::ClosedBall) = false
isempty(d::OpenBall) = radius(d) == 0

==(d1::Ball, d2::Ball) = isclosedset(d1)==isclosedset(d2) &&
    radius(d1)==radius(d2) && center(d1)==center(d2)
hash(d::Ball, h::UInt) = hashrec("Ball", isclosedset(d), radius(d), center(d), h)

function issubset(d1::Ball, d2::Ball)
    if dimension(d1) == dimension(d2)
        if center(d1) == center(d2)
            if radius(d1) < radius(d2)
                true
            elseif radius(d1) == radius(d2)
                isclosedset(d2) || isopenset(d1)
            else
                false
            end
        else
            dist = norm(center(d1)-center(d2))
            if dist + radius(d1) < radius(d2)
                true
            else
                false
            end
        end
    else
        false
    end
end

normal(d::Ball, x) = normal(boundary(d), x)

distance_to(d::Ball, x) = x ∈ d ? zero(prectype(d)) : norm(x-center(d))-radius(d)

# We choose the center of the ball here. Concrete types should implement 'center'
point_in_domain(d::Ball) = center(d)


"The unit ball."
abstract type UnitBall{T,C} <: Ball{T,C} end

radius(::UnitBall) = 1
center(d::UnitBall{T}) where {T<:StaticTypes} = zero(T)

UnitBall(n::Int) = DynamicUnitBall(n)
UnitBall(::Val{N} = Val(3)) where {N} = EuclideanUnitBall{N}()

UnitBall{T}(n::Int) where {T <: StaticTypes} = StaticUnitBall{T}(n)
UnitBall{T}(::Val{N}) where {N,T} = StaticUnitBall{T}(Val(N))
UnitBall{T}() where {T <: StaticTypes} = StaticUnitBall{T}()
UnitBall{T}(n::Int) where {T} = DynamicUnitBall{T}(n)

UnitBall{T,C}(n::Int) where {T <: StaticTypes,C} = StaticUnitBall{T,C}(n)
UnitBall{T,C}(::Val{N}) where {N,T,C} = StaticUnitBall{T,C}(Val(N))
UnitBall{T,C}() where {T <: StaticTypes,C} = StaticUnitBall{T,C}()
UnitBall{T,C}(n::Int) where {T,C} = DynamicUnitBall{T,C}(n)

const ClosedUnitBall{T} = UnitBall{T,:closed}
const OpenUnitBall{T} = UnitBall{T,:open}

indomain(x, d::OpenUnitBall) = norm(x) < 1
indomain(x, d::ClosedUnitBall) = norm(x) <= 1
approx_indomain(x, d::OpenUnitBall, tolerance) = norm(x) < 1+tolerance
approx_indomain(x, d::ClosedUnitBall, tolerance) = norm(x) <= 1+tolerance



==(d1::UnitBall, d2::UnitBall) = isclosedset(d1)==isclosedset(d2) &&
    dimension(d1) == dimension(d2)

issubset(d1::UnitBall, d2::UnitBall) =
    dimension(d1) == dimension(d2) && (isclosedset(d2) || isopenset(d1))

convert(::Type{SublevelSet}, d::UnitBall{T,C}) where {T,C} =
    SublevelSet{T,C}(norm, radius(d))
convert(::Type{SublevelSet{T}}, d::UnitBall{S,C}) where {T,S,C} =
    SublevelSet{T,C}(norm, radius(d))


"The unit ball with fixed dimension(s) specified by the element type."
struct StaticUnitBall{T,C} <: UnitBall{T,C}
end

StaticUnitBall() = StaticUnitBall{Float64}()
StaticUnitBall(::Val{N}) where {N} = StaticUnitBall{SVector{N,Float64}}()

StaticUnitBall{T}() where {T} = StaticUnitBall{T,:closed}()

StaticUnitBall{T}(n::Int) where {T} =
    (@assert n == euclideandimension(T); StaticUnitBall{T}())
StaticUnitBall{T}(::Val{N}) where {N,T} =
    (@assert N == euclideandimension(T); StaticUnitBall{T}())

StaticUnitBall{T,C}(n::Int) where {T,C} =
    (@assert n == euclideandimension(T); StaticUnitBall{T,C}())
StaticUnitBall{T,C}(::Val{N}) where {N,T,C} =
    (@assert N == euclideandimension(T); StaticUnitBall{T,C}())


"The unit ball in a fixed N-dimensional space."
const EuclideanUnitBall{N,T,C} = StaticUnitBall{SVector{N,T},C}

EuclideanUnitBall{N}() where {N} = EuclideanUnitBall{N,Float64}()

similardomain(d::StaticUnitBall{S,C}, ::Type{T}) where {S,C,T<:StaticTypes} =
    StaticUnitBall{T,C}()
similardomain(d::StaticUnitBall{S,C}, ::Type{T}) where {S,C,T} =
    DynamicUnitBall{T,C}(dimension(d))

const UnitDisk{T} = UnitBall{SVector{2,T},:closed}
UnitDisk() = UnitDisk{Float64}()

## A StaticUnitBall{<:Number} equals the interval [-1,1]  (open or closed)
convert(::Type{Interval}, d::StaticUnitBall{T,:closed}) where {T <: Number} =
    ChebyshevInterval{T}()
convert(::Type{Interval}, d::StaticUnitBall{T,:open}) where {T <: Number} =
    OpenInterval{T}(-1, 1)

canonicaldomain(::Equal, d::StaticUnitBall{T}) where {T<:Number} = convert(Interval, d)


"The unit ball with variable dimension stored in a data field."
struct DynamicUnitBall{T,C} <: UnitBall{T,C}
    dimension   ::  Int

    DynamicUnitBall{T,C}(n::Int) where {T,C} = new(n)
    DynamicUnitBall{T,C}(n::Int) where {T<:StaticTypes,C} =
        (@assert n == euclideandimension(T); new(n))
end

DynamicUnitBall(n::Int) = DynamicUnitBall{Vector{Float64}}(n)
DynamicUnitBall{T}(n::Int) where {T} = DynamicUnitBall{T,:closed}(n)

dimension(d::DynamicUnitBall) = d.dimension

center(d::DynamicUnitBall{T}) where {T} = zeros(eltype(T), dimension(d))

"The unit ball with vector elements of a given dimension."
const VectorUnitBall{T,C} = DynamicUnitBall{Vector{T},C}

VectorUnitBall(n::Int = 3) = VectorUnitBall{Float64}(n)
VectorUnitDisk() = VectorUnitBall(2)

similardomain(d::DynamicUnitBall{S,C}, ::Type{T}) where {S,C,T} =
    DynamicUnitBall{T,C}(dimension(d))
similardomain(d::DynamicUnitBall{S,C}, ::Type{T}) where {S,C,T<:StaticTypes} =
    StaticUnitBall{T,C}()


"A `GenericBall` is a ball with a given radius and center."
struct GenericBall{T,C,S} <: Ball{T,C}
    radius  ::  S
    center  ::  T
end

GenericBall() = GenericBall(1.0)
GenericBall(radius::S) where {S<:Number} = GenericBall(radius, zero(SVector{3,float(S)}))
GenericBall(radius::S, center::T) where {S<:Number,U,T<:Vector{U}} =
    GenericBall{Vector{promote_type(S,U)}}(radius, center)
GenericBall(radius::S, center::T) where {S<:Number,N,U,T<:SVector{N,U}} =
    GenericBall{SVector{N,promote_type(S,U)}}(radius, center)
GenericBall(radius::S, center::T) where {S<:Number,T<:AbstractVector} =
    GenericBall{Vector{promote_type(S,eltype(T))}}(radius, center)
GenericBall(radius::S, center::T) where {S<:Number,T<:Number} =
    GenericBall{promote_type(S,T)}(radius, center)
GenericBall{T}(radius) where {T} = GenericBall{T,:closed}(radius)
GenericBall{T}(radius, center) where {T} = GenericBall{T,:closed}(radius, center)
GenericBall{T,C}(radius) where {T,C} = GenericBall{T,C}(radius, zero(T))
GenericBall{T,C}(radius, center) where {T,C} = GenericBall{T,C,eltype(T)}(radius, center)

radius(d::GenericBall) = d.radius
center(d::GenericBall) = d.center

dimension(d::GenericBall) = length(d.center)

canonicaldomain(d::GenericBall{T,C}) where {T,C} = UnitBall{T,C}(dimension(d))
mapfrom_canonical(d::GenericBall) = AffineMap(radius(d), center(d))

boundingbox(d::GenericBall) =
    map_boundingbox(boundingbox(canonicaldomain(d)), mapfrom_canonical(d))

# Preserve the `Ball` type under affine maps which preserve shape
map_domain(m::GenericAffineMap{T,S}, d::Ball{U,C}) where {T<:AbstractVector,S<:Number,U<:AbstractVector,C} =
    Ball{U,C}(m.A * radius(d), m.A * center(d) + m.b)
map_domain(m::AffineMap{T}, d::Ball{U,C}) where {T<:Number,U<:Number,C} =
    Ball{promote_type(U,T),C}(m.A * radius(d), m.A * center(d) + m.b)
map_domain(m::GenericLinearMap{T,S}, d::Ball{U,C}) where {S<:Number,T<:AbstractVector,U<:AbstractVector,C} =
    Ball{promote_type(T,U),C}(m.A * radius(d), m.A * center(d))
map_domain(m::LinearMap{T}, d::Ball{U,C}) where {T<:Number,U<:Number,C} =
    Ball{promote_type(T,U),C}(m.A * radius(d), m.A * center(d))
map_domain(m::LinearMap{T}, d::Ball{U,C}) where {T<:Number,U,C} =
    Ball{U,C}(m.A * radius(d), m.A * center(d))
map_domain(m::Translation{T}, d::Ball{S,C}) where {T<:AbstractVector,S<:AbstractVector,C} =
    Ball{S,C}(radius(d), center(d) + m.b)
map_domain(m::Translation{T}, d::Ball{S,C}) where {T<:Number,S<:Number,C} =
    Ball{S,C}(radius(d), center(d) + m.b)


show(io::IO, d::EuclideanUnitBall{3,Float64,:closed}) = print(io, "UnitBall()")
show(io::IO, d::EuclideanUnitBall{N,Float64,:closed}) where {N} = print(io, "UnitBall(Val($(N)))")
show(io::IO, d::EuclideanUnitBall{3,Float64,:open}) = print(io, "UnitBall()  (open)")
show(io::IO, d::EuclideanUnitBall{N,Float64,:open}) where {N} = print(io, "UnitBall(Val($(N)))  (open)")
Display.object_parentheses(d::EuclideanUnitBall{N,Float64,:open}) where {N} = true
show(io::IO, d::EuclideanUnitBall{2,Float64,:closed}) = print(io, "UnitDisk()")
show(io::IO, d::UnitDisk{T}) where {T} = print(io, "UnitDisk{$(T)}()")
show(io::IO, d::UnitBall{T}) where {T<:Number} = print(io, "UnitBall{$(T)}()")
show(io::IO, d::UnitBall{T,:open}) where {T<:Number} = print(io, "UnitBall{$(T)}()  (open)")
Display.object_parentheses(d::UnitBall{T,:open}) where {T<:Number} = true
show(io::IO, d::VectorUnitBall{Float64,:closed}) = print(io, "UnitBall($(dimension(d)))")
show(io::IO, d::VectorUnitBall{Float64,:open}) = print(io, "UnitBall($(dimension(d)))  (open)")
show(io::IO, d::Ball) = print(io, "Ball($(radius(d)), $(center(d)))")


# The type hierarchy of spheres parallels that of Ball above:
# abstract Sphere
# |-> abstract UnitSphere: radius is 1
#     |-> StaticUnitSphere: dimension is part of type
#     |-> DynamicUnitSphere: dimension is specified by int field
# |-> GenericSphere: stores center and radius. Here, the dimension can be
#     obtained from the center, so there is only one type.
# There are aliases for SVector{N,T} and Vector{T}.

"Supertype of spherical domains for which elements satisfy `norm(x) == radius(sphere)`."
abstract type Sphere{T} <: Domain{T} end

indomain(x, d::Sphere) = norm(x-center(d)) == radius(d)
approx_indomain(x, d::Sphere, tolerance) =
    radius(d)-tolerance <= norm(x-center(d)) <= radius(d)+tolerance

# Convenience constructors for the abstract type
Sphere() = UnitSphere()
Sphere{T}() where {T} = UnitSphere{T}()
Sphere(radius::Number) = GenericSphere(radius)
Sphere{T}(radius::Number) where {T} = GenericSphere{T}(radius)
Sphere(radius::Number, center) = GenericSphere(radius, center)
Sphere{T}(radius::Number, center) where {T} = GenericSphere{T}(radius, center)

isempty(::Sphere) = false

isclosedset(::Sphere) = true
isopenset(::Sphere) = false

normal(d::Sphere, x) = (x-center(x))/norm(x-center(x))

distance_to(d::Sphere, x) = abs(norm(x-center(d))-radius(d))

point_in_domain(d::Sphere) = center(d) + unitvector(d, 1)

==(d1::Sphere, d2::Sphere) = radius(d1)==radius(d2) && center(d1)==center(d2)

"A hypersphere in a fixed N-dimensional Euclidean space."
const EuclideanSphere{N,T} = Sphere{SVector{N,T}}

"A sphere with vector elements of variable length."
const VectorSphere{T} = Sphere{Vector{T}}


"The unit sphere."
abstract type UnitSphere{T} <: Sphere{T} end

radius(::UnitSphere) = 1
center(d::UnitSphere{T}) where {T<:StaticTypes} = zero(T)

normal(d::UnitSphere, x) = x/norm(x)

indomain(x, d::UnitSphere) = norm(x) == 1
approx_indomain(x, d::UnitSphere, tolerance) =
    1-tolerance <= norm(x) <= 1+tolerance


UnitSphere(n::Int) = DynamicUnitSphere(n)
UnitSphere(::Val{N} = Val(3)) where {N} = EuclideanUnitSphere{N}()

UnitSphere{T}(n::Int) where {T <: StaticTypes} = StaticUnitSphere{T}(n)
UnitSphere{T}(::Val{N}) where {N,T} = StaticUnitSphere{T}(Val(N))
UnitSphere{T}() where {T <: StaticTypes} = StaticUnitSphere{T}()
UnitSphere{T}(n::Int) where {T} = DynamicUnitSphere{T}(n)

==(d1::UnitSphere, d2::UnitSphere) = dimension(d1)==dimension(d2)

issubset1(d1::UnitSphere, d2::UnitBall) =
    dimension(d1) == dimension(d2) && isclosedset(d2)

convert(::Type{LevelSet}, d::UnitSphere{T}) where {T} = LevelSet{T}(norm, radius(d))
convert(::Type{LevelSet{T}}, d::UnitSphere) where {T} = LevelSet{T}(norm, radius(d))


"The unit sphere with fixed dimension(s) specified by the element type."
struct StaticUnitSphere{T} <: UnitSphere{T}
end

StaticUnitSphere() = StaticUnitSphere{Float64}()
StaticUnitSphere(::Val{N}) where {N} = StaticUnitSphere{SVector{N,Float64}}()

StaticUnitSphere{T}(n::Int) where {T} =
    (@assert n == euclideandimension(T); StaticUnitSphere{T}())
StaticUnitSphere{T}(::Val{N}) where {N,T} =
    (@assert N == euclideandimension(T); StaticUnitSphere{T}())

similardomain(d::StaticUnitSphere, ::Type{T}) where {T<:StaticTypes} =
    StaticUnitSphere{T}()
similardomain(d::StaticUnitSphere, ::Type{T}) where {T} =
    DynamicUnitSphere{T}(dimension(d))


"The unit sphere in a fixed N-dimensional Euclidean space."
const EuclideanUnitSphere{N,T} = StaticUnitSphere{SVector{N,T}}

EuclideanUnitSphere{N}() where {N} = EuclideanUnitSphere{N,Float64}()


"The unit circle in 2D."
const UnitCircle{T} = UnitSphere{SVector{2,T}}
UnitCircle() = UnitCircle{Float64}()

"The unit sphere with variable dimension."
struct DynamicUnitSphere{T} <: UnitSphere{T}
    dimension   ::  Int

    DynamicUnitSphere{T}(n::Int) where {T} = new(n)
    DynamicUnitSphere{T}(n::Int) where {T<:StaticTypes} =
        (@assert n == euclideandimension(T); new(n))
end

DynamicUnitSphere(n::Int) = DynamicUnitSphere{Vector{Float64}}(n)

dimension(d::DynamicUnitSphere) = d.dimension

center(d::DynamicUnitSphere{T}) where {T} = zeros(eltype(T), dimension(d))

"The unit sphere with vector elements of a given dimension."
const VectorUnitSphere{T} = DynamicUnitSphere{Vector{T}}

VectorUnitSphere(n::Int = 3) = VectorUnitSphere{Float64}(n)
VectorUnitCircle() = VectorUnitSphere(2)

similardomain(d::DynamicUnitSphere, ::Type{T}) where {T} =
    DynamicUnitSphere{T}(d.dimension)
similardomain(d::DynamicUnitSphere, ::Type{T}) where {T <: StaticTypes} =
    StaticUnitSphere{T}()



"A `GenericSphere` is a sphere with a given radius and center."
struct GenericSphere{T,S} <: Sphere{T}
    radius  ::  S
    center  ::  T
end

GenericSphere() = GenericSphere(1.0)
GenericSphere(radius::S) where {S<:Number} =
    GenericSphere(radius, zero(SVector{3,float(S)}))
GenericSphere(radius::S, center::T) where {S<:Number,U,T<:Vector{U}} =
    GenericSphere{Vector{promote_type(S,U)}}(radius, center)
GenericSphere(radius::S, center::T) where {S<:Number,N,U,T<:SVector{N,U}} =
    GenericSphere{SVector{N,promote_type(S,U)}}(radius, center)
GenericSphere(radius::S, center::T) where {S<:Number,T<:AbstractVector} =
    GenericSphere{Vector{promote_type(S,eltype(T))}}(radius, center)
GenericSphere{T}(radius) where {T} = GenericSphere{T}(radius, zero(T))
GenericSphere{T}(radius, center) where {T} =
    GenericSphere{T,eltype(T)}(radius, center)

radius(d::GenericSphere) = d.radius
center(d::GenericSphere) = d.center

dimension(d::GenericSphere) = length(d.center)

canonicaldomain(d::GenericSphere{T}) where {T} = UnitSphere{T}(dimension(d))
mapfrom_canonical(d::GenericSphere) = AffineMap(radius(d), center(d))

boundingbox(d::GenericSphere) =
    map_boundingbox(boundingbox(canonicaldomain(d)), mapfrom_canonical(d))

# Preserve the `Sphere` type under affine maps which preserve shape
map_domain(m::GenericAffineMap{T,S}, d::Sphere{U}) where {T<:AbstractVector,S<:Number,U<:AbstractVector} =
    Sphere{U}(m.A * radius(d), m.A * center(d) + m.b)
map_domain(m::AffineMap{T}, d::Sphere{U}) where {T<:Number,U<:Number} =
    Sphere{promote_type(U,T)}(m.A * radius(d), m.A * center(d) + m.b)
map_domain(m::GenericLinearMap{T,S}, d::Sphere{U}) where {S<:Number,T<:AbstractVector,U<:AbstractVector} =
    Sphere{promote_type(T,U)}(m.A * radius(d), m.A * center(d))
map_domain(m::LinearMap{T}, d::Sphere{U}) where {T<:Number,U<:Number} =
    Sphere{promote_type(T,U)}(m.A * radius(d), m.A * center(d))
map_domain(m::LinearMap{T}, d::Sphere{U}) where {T<:Number,U} =
    Sphere{U}(m.A * radius(d), m.A * center(d))
map_domain(m::Translation{T}, d::Sphere{S}) where {T<:AbstractVector,S<:AbstractVector} =
    Sphere{S}(radius(d), center(d) + m.b)
map_domain(m::Translation{T}, d::Sphere{S}) where {T<:Number,S<:Number} =
    Sphere{S}(radius(d), center(d) + m.b)

show(io::IO, d::EuclideanUnitSphere{3,Float64}) = print(io, "UnitSphere()")
show(io::IO, d::EuclideanUnitSphere{N,Float64}) where {N} = print(io, "UnitSphere(Val($(N)))")
show(io::IO, d::EuclideanUnitSphere{2,Float64}) = print(io, "UnitCircle()")
show(io::IO, d::UnitCircle{T}) where {T} = print(io, "UnitCircle{$(T)}()")
show(io::IO, d::VectorUnitSphere{Float64}) = print(io, "UnitSphere($(dimension(d)))")
show(io::IO, d::UnitSphere{T}) where {T<:Number} = print(io, "UnitSphere{$(T)}()")
show(io::IO, d::Sphere) = print(io, "Sphere($(radius(d)), $(center(d)))")


boundary(d::UnitBall{T}) where {T} = UnitSphere{T}(dimension(d))
boundary(d::Ball{T}) where {T} = Sphere{T}(radius(d), center(d))

interior(d::StaticUnitBall{T}) where {T} = StaticUnitBall{T,:open}()
interior(d::DynamicUnitBall{T}) where {T} = DynamicUnitBall{T,:open}(dimension(d))
interior(d::Ball{T}) where {T} = Ball{T,:open}(radius(d), center(d))

closure(d::StaticUnitBall{T}) where {T} = StaticUnitBall{T,:closed}()
closure(d::DynamicUnitBall{T}) where {T} = DynamicUnitBall{T,:closed}(dimension(d))
closure(d::Ball{T}) where {T} = Ball{T,:closed}(radius(d), center(d))

boundingbox(d::UnitBall{T}) where {T<:Number} = ChebyshevInterval{T}()
boundingbox(d::UnitBall{SVector{N,T}}) where {N,T} =
    ChebyshevProductDomain{N,T}()
boundingbox(d::UnitBall{T}) where {T} =
    Rectangle{T}(-ones(eltype(T), dimension(d)), ones(eltype(T), dimension(d)))

boundingbox(d::UnitSphere{T}) where {T<:Number} = ChebyshevInterval{T}()
boundingbox(d::UnitSphere{SVector{N,T}}) where {N,T} =
    ChebyshevProductDomain{N,T}()
boundingbox(d::UnitSphere{T}) where {T} =
    Rectangle{T}(-ones(eltype(T), dimension(d)), ones(eltype(T), dimension(d)))



################
# Applications
################

"Create a cylinder with given radius and length."
cylinder(::Type{T} = Float64) where {T} = UnitDisk{T}() × UnitInterval{T}()
cylinder(radius, length) = cylinder(promote(radius, length)...)
cylinder(radius::T, length::T) where {T} = (radius .* UnitDisk{T}()) × (0..length)

"Create an ellipse curve with semi-axes lengths `a` and `b` respectively."
ellipse(a::Number, b::Number) = ellipse(promote(a,b)...)
ellipse(a::T, b::T) where {T <: Number} = LinearMap(a, b).(UnitCircle{T}())

"Create an ellipse-shaped domain with semi-axes lengths `a` and `b` respectively."
ellipse_shape(a::Number, b::Number) = ellipse_shape(promote(a,b)...)
ellipse_shape(a::T, b::T) where {T <: Number} = LinearMap(a, b).(UnitDisk{T}())


"""
The map `[cos(2πt), sin(2πt)]` from `[0,1]` to the unit circle in `ℝ^2`.
"""
struct UnitCircleMap{T} <: Map{T} end

UnitCircleMap() = UnitCircleMap{Float64}()

mapsize(m::UnitCircleMap) = (2,)

applymap(m::UnitCircleMap{T}, t) where {T} = SVector(cos(2*T(pi)*t), sin(2*T(pi)*t))

function jacobian(m::UnitCircleMap{T}, t) where {T}
    a = 2*T(pi)
    SVector(-a*sin(a*t), a*cos(a*t))
end

# we know that the differential volume is 2*pi
diffvolume(m::UnitCircleMap{T}) where {T} = ConstantMap{T}(2*T(pi))

"""
`AngleMap` is a left inverse of `UnitCircleMap`. A 2D vector `x` is projected onto
the intersection point with the unit circle of the line connecting `x` to the
origin. The angle of this point, scaled to the interval `[0,1)`, is the result.
"""
struct AngleMap{T} <: Map{SVector{2,T}}
end

AngleMap() = AngleMap{Float64}()

function applymap(m::AngleMap{T}, x) where {T}
    twopi = 2*convert(T, pi)
    θ = atan(x[2],x[1])
    if θ < 0
        # atan2 returns an angle in (-π,π], convert to [0,2π) using periodicity.
        θ += twopi
    end
    # And divide by 2π to scale to [0,1)
    θ / twopi
end

mapsize(m::AngleMap) = (1,2)

function jacobian(m::AngleMap{T}, t) where {T}
    x = t[1]; y = t[2]
    twopi = 2*convert(T, pi)
    v = SVector(one(T) / (1+(y/x)^2) * (-y/x^2) * 1/twopi,
            one(T) / (1+(y/x)^2) * one(T)/x * 1/twopi)
    transpose(v)
end


leftinverse(m::UnitCircleMap{T}) where {T} = AngleMap{T}()
leftinverse(m::UnitCircleMap, x) = leftinverse(m)(x)
rightinverse(m::AngleMap{T}) where {T} = UnitCircleMap{T}()
rightinverse(m::AngleMap, x) = rightinverse(m)(x)

canonicaldomain(::Parameterization, d::UnitCircle{T}) where {T} = UnitInterval{T}()
mapfrom_canonical(::Parameterization, d::UnitCircle{T}) where {T} = UnitCircleMap{T}()


## The complex plane

const ComplexUnitCircle{T} = StaticUnitSphere{Complex{T}}
const ComplexUnitDisk{T,C} = StaticUnitBall{Complex{T},C}

ComplexUnitCircle() = ComplexUnitCircle{Float64}()
ComplexUnitDisk() = ComplexUnitDisk{Float64}()
ComplexUnitDisk{T}() where {T} = ComplexUnitDisk{T,:closed}()

show(io::IO, d::ComplexUnitCircle{Float64}) = print(io, "ComplexUnitCircle()")
show(io::IO, d::ComplexUnitDisk{Float64,:closed}) = print(io, "ComplexUnitDisk()")
show(io::IO, d::ComplexUnitDisk{Float64,:open}) = print(io, "ComplexUnitDisk()  (open)")
show(io::IO, d::ComplexUnitDisk{T,:closed}) where {T} = print(io, "ComplexUnitDisk{$(T)}()")
show(io::IO, d::ComplexUnitDisk{T,:open}) where {T} = print(io, "ComplexUnitDisk{$(T)}()  (open)")

canonicaldomain(::Parameterization, d::ComplexUnitCircle{T}) where {T} =
    UnitInterval{T}()
mapfrom_canonical(::Parameterization, d::ComplexUnitCircle{T}) where {T} =
    VectorToComplex{T}() ∘ UnitCircleMap{T}()
