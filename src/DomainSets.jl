module DomainSets

using StaticArrays
using LinearAlgebra, Statistics
import LinearAlgebra: cross, ×, pinv

using IntervalSets

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
    size, length, ndims, getindex, eltype, ndims, hash,
    inv, isreal, zero,
    # Types, promotions and conversions
    convert, widen, promote_rule,
    # Display
    show,
    # Various
    Bool, angle

# IntervalSets
import IntervalSets: (..), endpoints, Domain, AbstractInterval, TypedEndpointsInterval,
                        leftendpoint, rightendpoint, isleftopen, isleftclosed,
                        isrightopen, isrightclosed, isopenset, isclosedset,
                        infimum, supremum
export ..


################################
## Exhaustive list of exports
################################

## Utils

# from util/common.jl
export prectype, numtype,
    convert_numtype, promote_numtype,
    convert_prectype, promote_prectype,
    iscomposite, elements, element, numelements
# from util/products.jl
export cartesianproduct


## Maps

# from maps/map.jl
export AbstractMap, Map, TypedMap,
    applymap,
    domaintype, codomaintype,
    inverse, leftinverse, rightinverse,
    jacobian, jacdet
# from maps/lazy.jl
export WrappedMap
# from maps/composite.jl
export Composition, ∘
# from maps/product.jl
export ProductMap, tensorproduct
# from maps/basic.jl
export IdentityMap, VectorIdentityMap, ZeroMap, UnityMap, ConstantMap,
    isconstant, constant
# from maps/affine.jl
export AffineMap, Translation, LinearMap,
    matrix, vector,
    islinear, isaffine,
    interval_map, scaling_map


## Generic domains

# from generic/domain.jl
export Domain, EuclideanDomain, VectorDomain,
    dimension,
    approx_in,
    isopenset, isclosedset, iscompact,
    boundary, ∂,
    point_in_domain

# from generic/lazy.jl
export DerivedDomain, superdomain, WrappedDomain

# from generic/productdomain.jl
export ProductDomain, VcatDomain, VectorProductDomain, TupleProductDomain

# from generic/mapped_domain.jl
export MappedDomain,
    map_domain,
    mapped_domain,
    forward_map,
    inverse_map

# from generic/setoperations.jl
export UnionDomain, IntersectionDomain, DifferenceDomain

# from applications/rotations.jl
export rotate

export infimum, supremum


## Specific domains

# from domains/trivial.jl
export EmptySpace, FullSpace,
    ℤ, ℚ, ℝ, ℂ, ℝ1, ℝ2, ℝ3, ℝ4
# from domains/interval.jl
export AbstractInterval, Interval, UnitInterval, ChebyshevInterval,
    OpenInterval, ClosedInterval,
    leftendpoint, rightendpoint, isleftopen, isrightopen,
    cardinality,
    HalfLine, NegativeHalfLine
# from domains/simplex.jl
export EuclideanUnitSimplex, VectorUnitSimplex, UnitSimplex,
    center, corners
# from domains/point.jl
export Point
# from domains/ball.jl
export UnitCircle, VectorUnitCircle,
    UnitSphere, VectorUnitSphere,
    UnitDisk, VectorUnitDisk,
    UnitBall, VectorUnitBall,
    UnitHyperBall,  UnitHyperSphere,
    EuclideanUnitBall, EuclideanUnitSphere,
    ComplexUnitCircle, ComplexUnitDisk,
    ellipse, ellipse_shape, cylinder,
    parameterization, gradient
# from domains/segment.jl
export Segment, arclength

## Applications
# from applications/rotation.jl
export rotation_map,
    CartToPolarMap, PolarToCartMap

include("util/common.jl")
include("util/products.jl")

include("maps/map.jl")
include("maps/lazy.jl")
include("maps/composite.jl")
include("maps/product.jl")
include("maps/basic.jl")
include("maps/affine.jl")
include("maps/arithmetics.jl")

include("generic/domain.jl")
include("generic/lazy.jl")
include("generic/productdomain.jl")
include("generic/setoperations.jl")
include("generic/mapped_domain.jl")
include("generic/promotion.jl")
include("generic/arithmetics.jl")
include("generic/broadcast.jl")

include("domains/trivial.jl")
include("domains/point.jl")
include("domains/interval.jl")
include("domains/simplex.jl")
include("domains/ball.jl")
include("domains/segment.jl")

include("applications/rotation.jl")

end # module
