module DomainSets

using StaticArrays
using LinearAlgebra, Statistics
import LinearAlgebra: cross, ×, pinv

using IntervalSets
using CompositeTypes, CompositeTypes.Display, CompositeTypes.Indexing

# deprecations in v0.5
@deprecate IntersectionDomain IntersectDomain
@deprecate DifferenceDomain SetdiffDomain
@deprecate FlexibleUnitCube DynamicUnitCube
@deprecate FlexibleUnitSphere DynamicUnitSphere
@deprecate FlexibleUnitBall DynamicUnitBall
@deprecate Composition ComposedMap

@deprecate element component
@deprecate elements components
@deprecate numelements ncomponents


################################
## Exhaustive list of imports
################################

# Generated functions
import Base: @ncall

# Operator symbols
import Base: *, +, -, /, \, ^,
    |, &,
    ∪, ∩,
    ==, isapprox,
    ∘,
    # Set operations
    intersect, union, setdiff, in, isempty, minimum, maximum,
    issubset,
    # Arrays
    eltype, hash, isreal,
    # Types, promotions and conversions
    convert, promote,
    # Display
    show

# IntervalSets
import IntervalSets: (..), endpoints, Domain, AbstractInterval, TypedEndpointsInterval,
                        leftendpoint, rightendpoint, isleftopen, isleftclosed,
                        isrightopen, isrightclosed, isopenset, isclosedset,
                        infimum, supremum
export ..


import CompositeTypes: component, components


################################
## Exhaustive list of exports
################################

## Utils

# from util/common.jl
export prectype, numtype,
    convert_numtype, promote_numtype,
    convert_prectype, promote_prectype,
    iscomposite, component, components, ncomponents

## Maps

# from maps/map.jl
export AbstractMap, Map, TypedMap,
    applymap,
    domaintype, codomaintype,
    inverse, leftinverse, rightinverse,
    mapsize, jacobian, jacdet, diffvolume
# from maps/composite.jl
export ComposedMap, composedmap, ∘
# from maps/product.jl
export ProductMap, productmap
# from maps/basic.jl
export IdentityMap,
    StaticIdentityMap, VectorIdentityMap,
    ZeroMap, UnityMap, ConstantMap,
    isconstant, constant
# from maps/affine.jl
export AffineMap, Translation, LinearMap,
    matrix, vector,
    islinear, isaffine


## Generic domains

# from generic/domain.jl
export Domain, EuclideanDomain, VectorDomain,
    dimension,
    approx_in,
    isopenset, isclosedset, iscompact,
    boundary, ∂,
    boundingbox,
    interior, closure,
    volume,
    point_in_domain,
    normal, tangents, distance_to,
    canonicaldomain, mapto_canonical, mapfrom_canonical, hascanonicaldomain,
    mapto,
    parameterdomain, parameterization, hasparameterization,
    mapfrom_parameterdomain, mapto_parameterdomain

# from generic/lazy.jl
export superdomain

# from generic/productdomain.jl
export ProductDomain, productdomain,
    VcatDomain, VectorProductDomain, TupleProductDomain

# from generic/mapped.jl
export MappedDomain,
    map_domain,
    mapped_domain,
    forward_map,
    inverse_map

# from generic/setoperations.jl
export UnionDomain, uniondomain,
    IntersectDomain, intersectdomain,
    SetdiffDomain, setdiffdomain

# from applications/rotations.jl
export rotate

export infimum, supremum


## Specific domains

# from domains/trivial.jl
export EmptySpace, FullSpace,
    emptyspace, fullspace,
    ℤ, ℚ, ℝ, ℂ, ℝ1, ℝ2, ℝ3, ℝ4
# from domains/interval.jl
export AbstractInterval, Interval, UnitInterval, ChebyshevInterval,
    OpenInterval, ClosedInterval,
    leftendpoint, rightendpoint, isleftopen, isrightopen,
    cardinality,
    HalfLine, NegativeHalfLine
# from domains/simplex.jl
export UnitSimplex,
    StaticUnitSimplex, DynamicUnitSimplex,
    EuclideanUnitSimplex, VectorUnitSimplex,
    corners
# from domains/point.jl
export Point
# from domains/ball.jl
export Ball, UnitBall,
    center, radius,
    StaticUnitBall, DynamicUnitBall,
    EuclideanUnitBall, VectorUnitBall,
    Sphere, UnitSphere,
    StaticUnitSphere, DynamicUnitSphere,
    VectorUnitSphere, EuclideanUnitSphere,
    UnitDisk, VectorUnitDisk,
    UnitCircle, VectorUnitCircle,
    ComplexUnitCircle, ComplexUnitDisk,
    ellipse, ellipse_shape, cylinder
# from domains/cube.jl
export UnitCube,
    StaticUnitCube, DynamicUnitCube,
    EuclideanUnitCube, VectorUnitCube,
    UnitSquare, UnitCube,
    Rectangle
# from domain/levelset.jl
export LevelSet, ZeroSet,
    SublevelSet, SubzeroSet,
    SuperlevelSet, SuperzeroSet,
    pseudolevel
# from domain/indicator.jl
export IndicatorFunction

## Applications
# from applications/rotation.jl
export rotation_map,
    CartToPolarMap, PolarToCartMap

include("util/common.jl")

include("maps/map.jl")
include("maps/lazy.jl")
include("maps/inverse.jl")
include("maps/jacobian.jl")
include("maps/composite.jl")
include("maps/product.jl")
include("maps/isomorphism.jl")
include("maps/basic.jl")
include("maps/affine.jl")
include("maps/arithmetics.jl")

include("generic/domain.jl")
include("generic/canonical.jl")
include("generic/lazy.jl")
include("generic/productdomain.jl")
include("generic/setoperations.jl")
include("generic/mapped.jl")
include("generic/broadcast.jl")

include("domains/trivial.jl")
include("domains/levelset.jl")
include("domains/point.jl")
include("domains/interval.jl")
include("domains/simplex.jl")
include("domains/ball.jl")
include("domains/cube.jl")
include("domains/indicator.jl")
include("domains/boundingbox.jl")

include("applications/coordinates.jl")
include("applications/rotation.jl")

end # module
