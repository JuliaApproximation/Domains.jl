using StaticArrays, DomainSets, Test

using DomainSets:
    MappedDomain,
    similar_interval,
    GenericBall, GenericSphere

struct Basis3Vector <: StaticVector{3,Float64} end

Base.getindex(::Basis3Vector, k::Int) = k == 1 ? 1.0 : 0.0

const io = IOBuffer()
const textmime = MIME"text/plain"()

struct NamedBall <: DomainSets.DerivedDomain{SVector{2,Float64}}
    domain  ::  Domain{SVector{2,Float64}}

    NamedBall() = new(2UnitDisk())
end

@testset "specific domains" begin
    @testset "empty space" begin
        d1 = EmptySpace()
        show(io,d1)
        @test isempty(d1)
        @test String(take!(io)) == "{} (empty domain)"
        @test eltype(d1) == Float64
        @test convert(Domain{BigFloat}, d1) === EmptySpace{BigFloat}()
        @test 0.5 ∉ d1
        @test !approx_in(0.5, d1)
        @test d1 ∩ d1 == d1
        @test d1 ∪ d1 == d1
        @test d1 \ d1 == d1
        @test boundary(d1) == d1
        @test dimension(d1) == 1
        @test isclosedset(d1)
        @test isopenset(d1)
        @test interior(d1) == d1
        @test closure(d1) == d1
        @test boundingbox(d1) == d1
        @test d1 == 2..1
        @test 2..1 == d1
        d2 = 0..1
        @test d1 ∩ d2 == d1
        @test d2 ∩ d1 == d1
        @test d1 ∪ d2 == d2
        @test d2 ∪ d1 == d2
        @test d1 \ d2 == d1
        @test d2 \ d1 == d2
        @test d2 \ d2 == d1
        # Test some promotions
        @test EmptySpace{Float64}() ∪ (0..1) isa AbstractInterval{Float64}
        @test EmptySpace{Int}() ∩ (0..1.0) isa EmptySpace{Float64}
        @test EmptySpace{Int}() \ (0..1.0) isa EmptySpace{Float64}
        @test (0..1) \ EmptySpace{Float64}() isa AbstractInterval{Float64}

        d2 = EmptySpace(SVector{2,Float64})
        @test isempty(d2)
        @test SA[0.1,0.2] ∉ d2
        @test [0.1,0.2] ∉ d2
        @test 1:2 ∉ d2
        @test !approx_in(SA[0.1,0.2], d2)
        @test boundary(d2) == d2
        @test dimension(d2) == 2

        @test emptyspace(0..1) == EmptySpace{Int}()
        @test emptyspace([1,2]) == EmptySpace{Int}()

        m = LinearMap(2)
        @test map_domain(m, emptyspace(Int)) == EmptySpace{Int}()
        @test mapped_domain(m, emptyspace(Int)) == EmptySpace{Int}()
    end

    @testset "full space" begin
        d1 = FullSpace()
        @test d1 == FullSpace{Float64}()
        show(io,d1)
        @test String(take!(io)) == "{x} (full space)"
        @test convert(Domain{BigFloat}, d1) === FullSpace{BigFloat}()
        @test DomainSets.euclideanspace(Val{2}()) == FullSpace{SVector{2,Float64}}()
        @test 0.5 ∈ d1
        @test point_in_domain(d1) == 0
        @test d1 ∪ d1 == d1
        @test d1 ∩ d1 == d1
        @test isempty(d1) == false
        @test boundary(d1) == EmptySpace{Float64}()
        @test isclosedset(d1)
        @test isopenset(d1)
        @test interior(d1) == d1
        @test closure(d1) == d1
        @test dimension(d1) == 1
        @test boundingbox(d1) == d1
        @test DomainSets.isequal1(d1, d1)
        @test DomainSets.isequal2(d1, d1)
        d2 = 0..1
        @test d1 ∪ d2 == d1
        @test d1 ∩ d2 == d2
        @test d2 ∩ d1 == d2
        @test (0..1.0) \ FullSpace{Int}() isa EmptySpace{Float64}
        @test typeof(FullSpace(0..1) .+ 1) <: FullSpace
        @test typeof(FullSpace(0..1) * 3) <: FullSpace
        @test infimum(d1) == typemin(Float64)
        @test supremum(d1) == typemax(Float64)
        @test FullSpace{Int}() == FullSpace{Float64}()

        d2 = FullSpace{SVector{2,Float64}}()
        @test SA[0.1,0.2] ∈ d2
        @test approx_in(SA[0.1,0.2], d2)
        @test !isempty(d2)
        @test boundary(d2) == EmptySpace{SVector{2,Float64}}()

        @test d2 == Domain(SVector{2,Float64})
        @test d2 == convert(Domain,SVector{2,Float64})
        @test d2 == convert(Domain{SVector{2,Float64}}, SVector{2,Float64})

        @test fullspace(0..1) == FullSpace{Int}()
        @test fullspace([1,2]) == FullSpace{Int}()

        @test uniondomain(UnitDisk(), FullSpace{SVector{2,Float64}}()) == FullSpace{SVector{2,Float64}}()
    end

    @testset "point" begin
        d = Domain(1.0)
        @test d isa Point
        @test 1 ∈ d
        @test 1.1 ∉ d
        @test approx_in(1.1, d, 0.2)
        @test !approx_in(1.2, d, 0.1)
        @test !isempty(d)
        @test boundary(d) == d
        @test boundingbox(d) == 1.0..1.0
        @test infimum(d) == d.x
        @test supremum(d) == d.x
        @test isclosedset(d)
        @test !isopenset(d)
        @test dimension(d) == 1
        @test isempty(interior(d))
        @test closure(d) == d
        @test canonicaldomain(d) == Point(0.0)
        @test mapfrom_canonical(d) == Translation(1.0)

        @test distance_to(d, 0.5) == abs(0.5-d.x)

        @test d .+ 1 == Domain(2.0)
        @test 1 .+ d == Domain(2.0)
        @test 1 .- d == Domain(0.0)
        @test d .- 1 == Domain(0.0)
        @test 2d  == Domain(2.0)
        @test d * 2 == Domain(2.0)
        @test d / 2 == Domain(0.5)
        @test 2 \ d == Domain(0.5)

        d1 = Domain(Set([1,2,3]))
        d2 = Point(1) ∪ Point(2) ∪ Point(3)

        @test d1 == d2

        @test convert(Domain{Float64}, Point(1)) ≡ Point(1.0)
        @test Number(Point(1)) ≡ convert(Number, Point(1)) ≡ convert(Int, Point(1)) ≡ 1
        @test convert(Domain{Float64}, 1) isa Point{Float64}

        @test point_in_domain(Point(1)) == 1

        @test Point(1) + Point(2) == Point(3)
        @test Point(1) - Point(2) == Point(-1)

        @test 0.5 ∉ (0..1)\Point(0.5)
        @test (0..1) \ Point(0.5) isa  UnionDomain{Float64}
        @test (0..1) \ Point(0.0) == Interval{:open,:closed,Float64}(0,1)
        @test (0..1) \ Point(1.0) == Interval{:closed,:open,Float64}(0,1)
        @test (0..1) \ Point(2.0) == Interval{:closed,:closed,Float64}(0,1)
        @test (0..1) \ 2.0 == (0..1) \ Point(2.0)
        @test issubset(Point(1), (0..2))
        @test Point(0.5) \ (0..1) == EmptySpace{Float64}()
        @test Point(0.5) \ (1..2) == Point(0.5)

        pv = Point([1,2,3])
        @test dimension(pv)==3
        @test canonicaldomain(pv) == Point([0,0,0])
        @test mapfrom_canonical(pv) == Translation(pv.x)
    end

    @testset "intervals" begin
        T = Float64
        @testset "ClosedInterval{$T}" begin
            d = zero(T)..one(T)
            @test approx_in(-0.1, d, 0.2)
            @test approx_in(1.1, d, 0.2)
            @test !approx_in(-0.2, d, 0.1)
            @test !approx_in(1.2, d, 0.1)
            @test isclosedset(d)
            @test !isopenset(d)
            @test similar_interval(d, big(0), big(1)) == Interval{:closed,:closed,BigInt}(0,1)
            @test closure(d) == d
            @test isopenset(interior(d))

            @test iscompact(d)
            @test typeof(similar_interval(d, one(T), 2*one(T))) == typeof(d)

            @test leftendpoint(d) ∈ ∂(d)
            @test rightendpoint(d) ∈ ∂(d)

            @test boundary(d) == uniondomain(Point(zero(T)), Point(one(T)))
            @test corners(d) == [0,1]
            @test boundingbox(d) == d

            @test similar_interval(0..1, 0, big(1.0)) isa ClosedInterval{BigFloat}

            @test canonicaldomain(d) isa ChebyshevInterval{Float64}
            @test canonicaldomain(0..1) isa ChebyshevInterval{Float64}
        end
        @testset "UnitInterval{$T}" begin
            d = UnitInterval{T}()
            @test leftendpoint(d) == zero(T)
            @test rightendpoint(d) == one(T)
            @test minimum(d) == infimum(d) == leftendpoint(d)
            @test maximum(d) == supremum(d) == rightendpoint(d)

            @test d ∩ d === d
            @test d ∪ d === d
            @test d \ d === EmptySpace{T}()

            @test isclosedset(d)
            @test !isopenset(d)
            @test iscompact(d)

            @test convert(Domain, d) ≡ d
            @test convert(Domain{T}, d) ≡ d
            @test convert(AbstractInterval, d) ≡ d
            @test convert(AbstractInterval{T}, d) ≡ d
            @test convert(UnitInterval, d) ≡ d
            @test convert(UnitInterval{T}, d) ≡ d
            @test convert(Domain{Float64}, d) ≡ UnitInterval()
            @test convert(AbstractInterval{Float64}, d) ≡ UnitInterval()
            @test convert(UnitInterval{Float64}, d) ≡ UnitInterval()

            @test canonicaldomain(d) === d
        end
        @testset "ChebyshevInterval{$T}" begin
            d = ChebyshevInterval{T}()
            @test leftendpoint(d) == -one(T)
            @test rightendpoint(d) == one(T)
            @test minimum(d) == infimum(d) == leftendpoint(d)
            @test maximum(d) == supremum(d) == rightendpoint(d)

            # test for promotion
            @test intersect(0..1, ChebyshevInterval()) == 0..1

            @test d ∩ d === d
            @test d ∪ d === d
            @test d \ d === EmptySpace{T}()
            unit = UnitInterval{T}()
            @test d ∩ unit === unit
            @test unit ∩ d === unit
            @test d ∪ unit === d
            @test unit ∪ d === d
            @test unit \ d === EmptySpace{T}()

            @test isclosedset(d)
            @test !isopenset(d)
            @test iscompact(d)

            @test convert(Domain, d) ≡ d
            @test convert(Domain{T}, d) ≡ d
            @test convert(AbstractInterval, d) ≡ d
            @test convert(AbstractInterval{T}, d) ≡ d
            @test convert(ChebyshevInterval, d) ≡ d
            @test convert(ChebyshevInterval{T}, d) ≡ d
            @test convert(Domain{Float64}, d) ≡ ChebyshevInterval()
            @test convert(AbstractInterval{Float64}, d) ≡ ChebyshevInterval()
            @test convert(ChebyshevInterval{Float64}, d) ≡ ChebyshevInterval()

            @test canonicaldomain(d) === d

            show(io, ChebyshevInterval())
            @test String(take!(io)) == "-1.0..1.0 (Chebyshev)"
        end
        @testset "HalfLine{$T}" begin
            d = HalfLine{T}()
            @test leftendpoint(d) == zero(T)
            @test rightendpoint(d) == T(Inf)
            @test minimum(d) == infimum(d) == leftendpoint(d)
            @test supremum(d) == rightendpoint(d)
            @test_throws ArgumentError maximum(d)

            @test d ∩ d === d
            @test d ∪ d === d
            @test d \ d == EmptySpace{T}()
            unit = UnitInterval{T}()
            cheb = ChebyshevInterval{T}()
            @test d ∩ unit === unit
            @test unit ∩ d === unit
            @test d ∩ cheb === unit
            @test cheb ∩ d === unit
            @test d ∪ unit === d
            @test unit ∪ d === d
            @test unit \ d === EmptySpace{T}()

            @test !isclosedset(d)
            @test !isopenset(d)
            @test !iscompact(d)
            @test 1. ∈ d
            @test -1. ∉ d
            @test approx_in(-0.1, d, 0.5)
            @test !approx_in(-0.5, d, 0.1)
            @test similar_interval(d, T(0), T(Inf)) == d

            @test 2d isa MappedDomain
            @test -2d isa MappedDomain

            @test boundary(d) == Point(0)
            @test leftendpoint(d) ∈ ∂(d)
            @test rightendpoint(d) ∉ ∂(d)
        end
        @testset "NegativeHalfLine{$T}" begin
            d = NegativeHalfLine{T}()
            @test leftendpoint(d) == -T(Inf)
            @test rightendpoint(d) == zero(T)
            @test infimum(d) == leftendpoint(d)
            @test supremum(d) == rightendpoint(d)
            @test_throws ArgumentError minimum(d)
            @test_throws ArgumentError maximum(d)

            @test d ∩ d === d
            @test d ∪ d === d
            @test d \ d == EmptySpace{T}()
            unit = UnitInterval{T}()
            cheb = ChebyshevInterval{T}()
            halfline = HalfLine{T}()
            @test unit ∩ d === EmptySpace{T}()
            @test d ∩ unit === EmptySpace{T}()
            @test d ∩ halfline === EmptySpace{T}()
            @test halfline ∩ d === EmptySpace{T}()
            @test d ∪ halfline === FullSpace{T}()
            @test halfline ∪ d === FullSpace{T}()
            @test unit \ d === unit
            @test cheb \ d === unit
            @test halfline \ d === halfline
            @test d \ unit === d
            @test d \ halfline === d

            @test !isclosedset(d)
            @test isopenset(d)
            @test !iscompact(d)
            @test -1. ∈ d
            @test 1. ∉ d
            @test approx_in(0.5, d, 1.)
            @test !approx_in(0.5, d, 0.4)
            @test similar_interval(d, T(-Inf), T(0)) == d

            @test 2d isa MappedDomain
            @test -2d isa MappedDomain

            @test boundary(d) == Point(0)
            @test leftendpoint(d) ∉ ∂(d)
            @test rightendpoint(d) ∈ ∂(d)
        end

        @testset "OpenInterval{$T}" begin
            d = OpenInterval(0,1)
            @test isopenset(d)
            @test closure(d) == UnitInterval()

            @test leftendpoint(d) ∈ ∂(d)
            @test rightendpoint(d) ∈ ∂(d)
        end

        @testset "Integer intervals" begin
            d = 0..1
            @test leftendpoint(d) ∈ ∂(d)
            @test rightendpoint(d) ∈ ∂(d)

            d = Interval{:open,:closed}(0,1)
            @test leftendpoint(d) ∈ ∂(d)
            @test rightendpoint(d) ∈ ∂(d)
            @test closure(d) == 0..1

            d = Interval{:closed,:open}(0,1)
            @test leftendpoint(d) ∈ ∂(d)
            @test rightendpoint(d) ∈ ∂(d)
            @test closure(d) == 0..1
        end

        @testset "approximate in for open and closed intervals" begin
            @test !approx_in(1.0, Interval{:open,:open}(0,1), 0)
            @test !approx_in(1.0, Interval{:closed,:open}(0,1), 0)
            @test !approx_in(0.0, Interval{:open,:open}(0,1), 0)
            @test !approx_in(0.0, Interval{:open,:closed}(0,1), 0)
            @test approx_in(0.0, Interval{:closed,:closed}(0,1), 0)
            @test approx_in(1.0, Interval{:closed,:closed}(0,1), 0)

            @test (0..1) ≈ (1e-16..1)
            @test (-1..0) ≈ (-1..1e-16)
            @test isapprox(HalfLine(), 0..Inf)
        end

        @testset "mapping between intervals" begin
            @test canonicaldomain(UnitInterval()) == UnitInterval()
            m = mapto(2..3, ChebyshevInterval())
            @test isaffine(m)
            @test m(2) ≈ -1
            @test m(3) ≈ 1
            m2 = mapto(4.0..6, 2..3)
            @test isaffine(m2)
            @test m2(4) ≈ 2
            @test m2(6) ≈ 3
        end

        @test DomainSets.isinterval(0..1)
        @test !DomainSets.isinterval(UnitBall())

        @test typeof(UnitInterval{Float64}(0.0..1.0)) <: UnitInterval
        @test typeof(ChebyshevInterval{Float64}(-1.0..1.0)) <: ChebyshevInterval

        ## Some mappings preserve the interval structure
        # Translation
        d = zero(T)..one(T)

        @test -Interval{:closed,:open}(2,3) isa Interval{:open,:closed}

        d2 = d .+ one(T)
        @test typeof(d2) == typeof(d)
        @test leftendpoint(d2) == one(T)
        @test rightendpoint(d2) == 2*one(T)

        d2 = one(T) .+ d
        @test typeof(d2) == typeof(d)
        @test leftendpoint(d2) == one(T)
        @test rightendpoint(d2) == 2*one(T)

        d2 = d .- one(T)
        @test typeof(d2) == typeof(d)
        @test leftendpoint(d2) == -one(T)
        @test rightendpoint(d2) == zero(T)

        d2 = -d
        @test typeof(d2) == typeof(d)
        @test leftendpoint(d2) == -one(T)
        @test rightendpoint(d2) == zero(T)

        d2 = one(T) .- d
        @test d2 == d

        # translation for UnitInterval
        # Does a shifted unit interval return an interval?
        d = UnitInterval{T}()
        d2 = d .+ one(T)
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == one(T)
        @test rightendpoint(d2) == 2*one(T)

        d2 = one(T) .+ d
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == one(T)
        @test rightendpoint(d2) == 2*one(T)

        d2 = d .- one(T)
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == -one(T)
        @test rightendpoint(d2) == zero(T)

        d2 = -d
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == -one(T)
        @test rightendpoint(d2) == zero(T)

        d2 = one(T) .- d
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == zero(T)
        @test rightendpoint(d2) == one(T)


        # translation for ChebyshevInterval
        d = ChebyshevInterval{T}()
        d2 = d .+ one(T)
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == zero(T)
        @test rightendpoint(d2) == 2*one(T)

        d2 = one(T) .+ d
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == zero(T)
        @test rightendpoint(d2) == 2*one(T)

        d2 = d .- one(T)
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == -2one(T)
        @test rightendpoint(d2) == zero(T)

        @test -d == d

        d2 = one(T) .- d
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == zero(T)
        @test rightendpoint(d2) == 2one(T)


        # Scaling
        d = zero(T)..one(T)
        d3 = T(2) * d
        @test typeof(d3) == typeof(d)
        @test leftendpoint(d3) == zero(T)
        @test rightendpoint(d3) == T(2)

        d3 = d * T(2)
        @test typeof(d3) == typeof(d)
        @test leftendpoint(d3) == zero(T)
        @test rightendpoint(d3) == T(2)

        d = zero(T)..one(T)
        d4 = d / T(2)
        @test typeof(d4) == typeof(d)
        @test leftendpoint(d4) == zero(T)
        @test rightendpoint(d4) == T(1)/T(2)

        d4 = T(2) \ d
        @test typeof(d4) == typeof(d)
        @test leftendpoint(d4) == zero(T)
        @test rightendpoint(d4) == T(1)/T(2)

        # Equality
        @test ChebyshevInterval() == ClosedInterval(-1.0,1.0) == ClosedInterval(-1,1)
        @test ChebyshevInterval() ≠ Interval{:closed,:open}(-1.0,1.0)
        @test ChebyshevInterval() ≈ ClosedInterval(-1.0,1.0) ≈ Interval{:closed,:open}(-1.0,1.0)

        # Union and intersection of intervals
        i1 = zero(T)..one(T)
        i2 = one(T)/3 .. one(T)/2
        i3 = one(T)/2 .. 2*one(T)
        i4 = T(2) .. T(3)
        # - union of completely overlapping intervals
        du1 = uniondomain(i1, i2)
        @test typeof(du1) <: AbstractInterval
        @test leftendpoint(du1) == leftendpoint(i1)
        @test rightendpoint(du1) == rightendpoint(i1)
        @test uniondomain(0..1, 0.1..1.5) isa AbstractInterval{Float64}

        # - intersection of completely overlapping intervals
        du2 = intersectdomain(i1, i2)
        @test typeof(du2) <: AbstractInterval
        @test leftendpoint(du2) == leftendpoint(i2)
        @test rightendpoint(du2) == rightendpoint(i2)

        # - union of partially overlapping intervals
        du3 = uniondomain(i1, i3)
        @test typeof(du3) <: AbstractInterval
        @test leftendpoint(du3) == leftendpoint(i1)
        @test rightendpoint(du3) == rightendpoint(i3)

        @test uniondomain(OpenInterval(0,1), 0..2) == 0..2
        @test uniondomain(OpenInterval(0,1), OpenInterval(0,2)) == OpenInterval(0,2)
        @test uniondomain(1..2, 0..1.5) == 0..2.0
        @test uniondomain(1..2.5, 0.8..1.5) == 0.8..2.5
        @test uniondomain(1..2.5, 0.8..2.5) == 0.8..2.5
        @test uniondomain(OpenInterval(1,2.5), OpenInterval(0.8,2.5)) == OpenInterval(0.8,2.5)

        # - intersection of partially overlapping intervals
        du4 = intersectdomain(i1, i3)
        @test typeof(du4) <: AbstractInterval
        @test leftendpoint(du4) == leftendpoint(i3)
        @test rightendpoint(du4) == rightendpoint(i1)

        # - union of non-overlapping intervals
        du5 = UnionDomain(i1) ∪ UnionDomain(i4)
        @test typeof(du5) <: UnionDomain

        # - intersection of non-overlapping intervals
        du6 = intersectdomain(i1, i4)
        @test isempty(du6)

        # - setdiff of intervals
        d1 = -2one(T).. 2one(T)
        @test d1 \ (3one(T) .. 4one(T)) == d1
        @test d1 \ (zero(T) .. one(T)) == UnionDomain((-2one(T)..zero(T))) ∪ UnionDomain((one(T).. 2one(T)))
        @test d1 \ (zero(T) .. 3one(T)) == (-2one(T) .. zero(T))
        @test d1 \ (-3one(T) .. zero(T)) == (zero(T) .. 2one(T))
        @test d1 \ (-4one(T) .. -3one(T)) == d1
        @test d1 \ (-4one(T) .. 4one(T)) == EmptySpace{T}()
        @test setdiffdomain(d1, (3one(T) .. 4one(T))) == d1
        @test setdiffdomain(d1, (zero(T) .. one(T))) == UnionDomain((-2one(T)..zero(T))) ∪ UnionDomain((one(T).. 2one(T)))
        @test setdiffdomain(d1, (zero(T) .. 3one(T))) == (-2one(T) .. zero(T))
        @test setdiffdomain(d1, (-3one(T) .. zero(T))) == (zero(T) .. 2one(T))
        @test setdiffdomain(d1, (-4one(T) .. -3one(T))) == d1
        @test setdiffdomain(d1, (-4one(T) .. 4one(T))) == EmptySpace{T}()

        # mixed types
        @test setdiffdomain(0..1, 0.0..0.5) == 0.5..1

        @test setdiffdomain(d1, -3) == d1
        @test setdiffdomain(d1, -2) == Interval{:open,:closed}(-2one(T),2one(T))
        @test setdiffdomain(d1, 2one(T)) == Interval{:closed,:open}(-2one(T),2one(T))
        @test setdiffdomain(d1, zero(T)) == UnionDomain(Interval{:closed,:open}(-2one(T),zero(T))) ∪ UnionDomain(Interval{:open,:closed}(zero(T),2one(T)))

        # - empty interval
        @test isempty(one(T)..zero(T))
        @test zero(T) ∉ (one(T)..zero(T))
        @test isempty(Interval{:open,:open}(zero(T),zero(T)))
        @test zero(T) ∉ Interval{:open,:open}(zero(T),zero(T))
        @test isempty(Interval{:open,:closed}(zero(T),zero(T)))
        @test zero(T) ∉ Interval{:open,:closed}(zero(T),zero(T))
        @test isempty(Interval{:closed,:open}(zero(T),zero(T)))
        @test zero(T) ∉ Interval{:closed,:open}(zero(T),zero(T))

        d = one(T) .. zero(T)
        @test_throws ArgumentError minimum(d)
        @test_throws ArgumentError maximum(d)
        @test_throws ArgumentError infimum(d)
        @test_throws ArgumentError supremum(d)

        # Subset relations of intervals
        @test issubset((zero(T)..one(T)), (zero(T).. 2*one(T)))
        @test issubset((zero(T)..one(T)), (zero(T).. one(T)))
        @test issubset(OpenInterval(zero(T),one(T)), zero(T) .. one(T))
        @test !issubset(zero(T) .. one(T), OpenInterval(zero(T), one(T)))
        @test issubset(UnitInterval{T}(), ChebyshevInterval{T}())

        # - convert
        d = zero(T).. one(T)
        @test d ≡ Interval(zero(T), one(T))
        @test d ≡ ClosedInterval(zero(T), one(T))

        @test convert(Domain, d) ≡ d
        @test convert(Domain{Float32}, d) ≡ (0f0 .. 1f0)
        @test convert(Domain{Float64}, d) ≡ (0.0 .. 1.0)
        @test convert(Domain, zero(T)..one(T)) ≡ d
        @test convert(Domain{T}, zero(T)..one(T)) ≡ d
        @test convert(AbstractInterval, zero(T)..one(T)) ≡ d
        @test convert(AbstractInterval{T}, zero(T)..one(T)) ≡ d
        @test convert(Interval, zero(T)..one(T)) ≡ d
        @test Interval(zero(T)..one(T)) ≡ d
        @test convert(ClosedInterval, zero(T)..one(T)) ≡ d
        @test ClosedInterval(zero(T)..one(T)) ≡ d
        @test convert(ClosedInterval{T}, zero(T)..one(T)) ≡ d
        @test ClosedInterval{T}(zero(T)..one(T)) ≡ d


        @testset "conversion from other types" begin
            @test convert(Domain{T}, 0..1) ≡ d
            @test convert(AbstractInterval{T}, 0..1) ≡ d
            @test convert(ClosedInterval{T}, 0..1) ≡ d
            @test ClosedInterval{T}(0..1) ≡ d
        end
    end

    @testset "balls" begin
        # Test UnitBall constructor
        @test UnitBall(2) isa VectorUnitBall{Float64}
        @test UnitBall(Val(2)) isa EuclideanUnitBall{2,Float64}
        @test UnitBall{Float64}() isa StaticUnitBall{Float64}
        @test UnitBall{Float64}(1) isa StaticUnitBall{Float64}
        @test_throws AssertionError UnitBall{Float64}(2)
        @test UnitBall{SVector{2,Float64}}(Val(2)) isa EuclideanUnitBall{2,Float64}
        @test_throws AssertionError UnitBall{SVector{2,Float64}}(Val(3))
        @test UnitBall{SVector{2,Float64}}(2) isa EuclideanUnitBall{2,Float64}
        @test_throws AssertionError UnitBall{SVector{2,Float64}}(3)
        @test UnitBall{SVector{2,Float64}}() isa EuclideanUnitBall{2,Float64}
        @test UnitBall{Vector{Float64}}(2) isa VectorUnitBall{Float64}
        @test_throws MethodError UnitBall{Vector{Float64}}()

        @test UnitBall{Float64,:open}() isa StaticUnitBall{Float64,:open}
        @test UnitBall{Float64,:closed}(1) isa StaticUnitBall{Float64,:closed}
        @test_throws AssertionError UnitBall{Float64,:closed}(2)
        @test UnitBall{SVector{2,Float64},:open}(Val(2)) isa EuclideanUnitBall{2,Float64,:open}
        @test_throws AssertionError UnitBall{SVector{2,Float64},:closed}(Val(3))
        @test UnitBall{SVector{2,Float64},:open}(2) isa EuclideanUnitBall{2,Float64,:open}
        @test_throws AssertionError UnitBall{SVector{2,Float64},:closed}(3)
        @test UnitBall{SVector{2,Float64},:open}() isa EuclideanUnitBall{2,Float64,:open}
        @test UnitBall{Vector{Float64},:closed}(2) isa VectorUnitBall{Float64,:closed}
        @test_throws MethodError UnitBall{Vector{Float64},:open}()

        @test DynamicUnitBall{Float64}(1) isa DynamicUnitBall{Float64}
        @test_throws AssertionError DynamicUnitBall{Float64}(2)

        # and then the Ball constructor
        @test Ball() isa UnitBall
        @test Ball(1.0, 2.0) isa DomainSets.GenericBall{Float64}
        @test Ball{Float64}() isa UnitBall{Float64}
        @test Ball{Float64,:open}() isa UnitBall{Float64,:open}
        @test Ball(1.0) == DomainSets.GenericBall{SVector{3,Float64},:closed}(1.0)
        @test Ball{BigFloat}(1) == DomainSets.GenericBall{BigFloat,:closed}(big(1))
        @test Ball{BigFloat,:open}(1) == DomainSets.GenericBall{BigFloat,:open}(big(1))
        @test_throws MethodError Ball{Vector{Float64}}(1.0)
        @test Ball(1.0, [1,2,3]) isa DomainSets.GenericBall{Vector{Float64},:closed,Float64}
        @test Ball{Vector{Float64}}(1.0, [1,2,3]) isa DomainSets.GenericBall{Vector{Float64},:closed,Float64}

        @test map_domain(AffineMap(2, [1,2,3]), UnitBall()) isa DomainSets.GenericBall
        @test map_domain(AffineMap(2, 3), UnitBall{Float64}()) isa DomainSets.GenericBall
        @test map_domain(LinearMap{SVector{3,Float64}}(2), UnitBall()) isa DomainSets.GenericBall
        @test map_domain(LinearMap(2), UnitBall{Float64}()) isa DomainSets.GenericBall
        @test map_domain(LinearMap(2), UnitBall()) isa DomainSets.GenericBall
        @test map_domain(Translation([1,2,3]), UnitBall()) isa DomainSets.GenericBall
        @test map_domain(Translation(1), UnitBall{Float64}()) isa DomainSets.GenericBall

        @test StaticUnitBall() isa StaticUnitBall{Float64}
        @test StaticUnitBall(Val(2)) isa StaticUnitBall{SVector{2,Float64}}

        @test GenericBall() == GenericBall(1.0)
        @test GenericBall(2.0, 1:5) isa GenericBall{Vector{Float64},:closed,Float64}
        @test GenericBall(2, 1.0:5.0) isa GenericBall{Vector{Float64},:closed,Float64}

        @test issubset(Ball(1.0, 2.0), Ball(1.0, 2.0))
        @test issubset(Ball(1.0, 2.0), Ball(2.0, 2.0))
        @test !issubset(Ball(2.0, 2.0), Ball(1.0, 2.0))
        @test issubset(Ball{Float64,:open}(2.0, 2.0), Ball(2.0, 2.0))
        @test !issubset(Ball(2.0, 2.0), Ball{Float64,:open}(2.0, 2.0))
        @test issubset(Ball(1.0, 0.0), Ball(2.5, 1.0))
        @test !issubset(Ball(1.0, 0.0), Ball(0.5, 1.0))
        @test !issubset(UnitBall(Val(2)), UnitBall(Val(3)))

        D = UnitDisk()
        @test SA[1.,0.] ∈ D
        @test SA[1.,1.] ∉ D
        @test approx_in(SA[1.0,0.0+1e-5], D, 1e-4)
        @test !isempty(D)
        @test isclosedset(D)
        @test !isopenset(D)
        D2 = convert(Domain{SVector{2,BigFloat}}, D)
        @test eltype(D2) == SVector{2,BigFloat}
        @test boundary(D) == UnitCircle()
        @test dimension(D) == 2
        @test boundingbox(D) == ProductDomain(ChebyshevInterval(), ChebyshevInterval())
        @test normal(UnitDisk(), [sqrt(2)/2, sqrt(2)/2]) ≈ [sqrt(2)/2, sqrt(2)/2]

        @test boundingbox(UnitBall{Float64}()) == ChebyshevInterval()

        @test convert(SublevelSet, UnitDisk()) isa SublevelSet{SVector{2,Float64},:closed}
        @test convert(SublevelSet{SVector{2,Float64}}, UnitDisk()) isa SublevelSet{SVector{2,Float64},:closed}
        @test convert(SublevelSet, EuclideanUnitBall{2,Float64,:open}()) isa SublevelSet{SVector{2,Float64},:open}

        @test convert(Interval, UnitBall{Float64}()) === ChebyshevInterval()
        @test convert(Interval, UnitBall{Float64,:open}()) === OpenInterval(-1.0, 1.0)
        @test UnitBall{Float64}() == ChebyshevInterval()

        @test convert(Domain{SVector{2,Float64}}, UnitBall(2)) isa StaticUnitBall
        @test convert(Domain{Vector{Float64}}, UnitBall(Val(2))) isa DynamicUnitBall

        @test repr(UnitBall()) == "UnitBall()"
        @test repr(UnitBall(Val(4))) == "UnitBall(Val(4))"
        @test repr(EuclideanUnitBall{3,Float64,:open}()) == "UnitBall()  (open)"
        @test repr(EuclideanUnitBall{4,Float64,:open}()) == "UnitBall(Val(4))  (open)"
        @test repr(UnitDisk()) == "UnitDisk()"
        @test repr(UnitDisk{BigFloat}()) == "UnitDisk{BigFloat}()"
        @test repr(UnitBall{Float64}()) == "UnitBall{Float64}()"
        @test repr(UnitBall{Float64,:open}()) == "UnitBall{Float64}()  (open)"
        @test repr(VectorUnitBall{Float64}(4)) == "UnitBall(4)"
        @test repr(VectorUnitBall{Float64,:open}(4)) == "UnitBall(4)  (open)"
        @test repr(Ball(1.0,2.0)) == "Ball(1.0, 2.0)"

        D = EuclideanUnitBall{2,Float64,:open}()
        @test !in(SA[1.0,0.0], D)
        @test in(SA[1.0-eps(Float64),0.0], D)
        @test approx_in(SA[1.1,0.0], D, 0.2)
        @test !approx_in(SA[1.1,0.0], D, 0.01)
        @test SA[0.2,0.2] ∈ D
        @test !isclosedset(D)
        @test isopenset(D)
        @test boundary(D) == UnitCircle()

        D = 2UnitDisk()
        @test D isa DomainSets.GenericBall
        @test SA[1.4, 1.4] ∈ D
        @test SA[1.5, 1.5] ∉ D
        @test typeof(1.2 * D)==typeof(D * 1.2)
        @test SA[1.5,1.5] ∈ 1.2 * D
        @test SA[1.5,1.5] ∈ D * 1.2
        @test !isempty(D)

        @test canonicaldomain(D) == UnitDisk()
        @test matrix(mapfrom_canonical(D)) == [2 0; 0 2]
        @test vector(mapfrom_canonical(D)) == [0; 0]
        @test parameterdomain(D) == canonicaldomain(D)
        @test mapfrom_parameterdomain(D) == mapfrom_canonical(D)
        @test boundingbox(D) == (-2.0..2.0)^2

        D = 2UnitDisk() .+ SA[1.0,1.0]
        @test D isa DomainSets.GenericBall
        @test boundary(D) isa DomainSets.GenericSphere
        @test SA[2.4, 2.4] ∈ D
        @test SA[3.5, 2.5] ∉ D
        @test !isempty(D)
        @test isopenset(interior(D))
        @test isclosedset(closure(D))
        @test canonicaldomain(D) isa UnitDisk

        B = UnitBall()
        @test SA[1.,0.0,0.] ∈ B
        @test SA[1.,0.1,0.] ∉ B
        @test !isempty(B)
        @test isclosedset(B)
        @test !isopenset(B)
        @test boundary(B) == UnitSphere()
        @test isopenset(interior(B))
        @test isclosedset(closure(B))

        B = 2UnitBall()
        @test D isa DomainSets.GenericBall
        @test SA[1.9,0.0,0.0] ∈ B
        @test SA[0,-1.9,0.0] ∈ B
        @test SA[0.0,0.0,-1.9] ∈ B
        @test SA[1.9,1.9,0.0] ∉ B
        @test !isempty(B)

        B = 2.0UnitBall() .+ SA[1.0,1.0,1.0]
        @test D isa DomainSets.GenericBall
        @test SA[2.9,1.0,1.0] ∈ B
        @test SA[1.0,-0.9,1.0] ∈ B
        @test SA[1.0,1.0,-0.9] ∈ B
        @test SA[2.9,2.9,1.0] ∉ B
        @test !isempty(B)

        C = VectorUnitBall(4)
        @test [1, 0, 0, 0] ∈ C
        @test [0.0,0.1,0.2,0.1] ∈ C
        @test SA[0.0,0.1,0.2,0.1] ∈ C
        @test_logs (:warn, "`in`: incompatible combination of vector with length 2 and domain 'UnitBall(4)' with dimension 4. Returning false.") [0.0,0.1] ∉ C
        @test [0.0,1.1,0.2,0.1] ∉ C
        @test !isempty(C)
        @test isclosedset(C)
        @test !isopenset(C)
        @test boundary(C) == VectorUnitSphere(4)
        cheb = ChebyshevInterval()
        @test boundingbox(C) == ProductDomain([cheb, cheb, cheb, cheb])
        @test dimension(C) == 4
        @test isopenset(interior(C))

        D = VectorUnitBall{Float64,:open}(4)
        @test !in([1, 0, 0, 0], D)
        @test in([1-eps(Float64), 0, 0, 0], D)
        @test approx_in([1.1, 0, 0, 0], D, 0.2)
        @test !approx_in([1.1, 0, 0, 0], D, 0.01)
        @test !isempty(D)
        @test approx_in(SA[1.01,0.0,0.0,0.0], D, 0.05)

        E = Ball{Float64,:open}(2.0)
        @test 0.5 ∈ E
        @test approx_in(0.5, E)
        @test approx_in(0.5, closure(E))
        @test isclosedset(StaticUnitBall{SVector{2,Float64}}())
        @test EuclideanUnitBall{2}() isa EuclideanUnitBall{2,Float64}

        show(io, EuclideanUnitBall{2,Float64,:open}())
        @test String(take!(io)) == "UnitBall(Val(2))  (open)"
        show(io, UnitCircle())
        @test String(take!(io)) == "UnitCircle()"
    end

    @testset "custom named ball" begin
        B = NamedBall()
        @test SA[1.4, 1.4] ∈ B
        @test SA[1.5, 1.5] ∉ B
        @test typeof(1.2 * B)==typeof(B * 1.2)
        @test SA[1.5,1.5] ∈ 1.2 * B
        @test SA[1.5,1.5] ∈ B * 1.2
        @test eltype(B) == eltype(2UnitDisk())
    end

    @testset "complex unit circle/disk" begin
        C = ComplexUnitCircle()
        @test eltype(C) == Complex{Float64}
        @test isclosedset(C)
        @test !isopenset(C)
        @test 1 ∈ C
        @test 1im ∈ C
        @test 1.1im ∉ C
        @test 0.2+0.5im ∉ C
        @test 1.2+0.5im ∉ C
        @test parameterdomain(C) == UnitInterval()

        D = ComplexUnitDisk()
        @test eltype(D) == Complex{Float64}
        @test isclosedset(D)
        @test !isopenset(D)
        @test 1 ∈ D
        @test 1im ∈ D
        @test 1.1im ∉ D
        @test 0.2+0.5im ∈ D
        @test 1.2+0.5im ∉ D

        D2 = ComplexUnitDisk{BigFloat,:open}()
        @test eltype(D2) == Complex{BigFloat}
        @test isopenset(D2)
        @test 1im ∉ D2
        @test 0.999 ∈ D2

        @test repr(ComplexUnitCircle()) == "ComplexUnitCircle()"
        @test repr(ComplexUnitDisk()) == "ComplexUnitDisk()"
        @test repr(ComplexUnitDisk{Float64,:open}()) == "ComplexUnitDisk()  (open)"
        @test repr(ComplexUnitDisk{BigFloat}()) == "ComplexUnitDisk{BigFloat}()"
        @test repr(ComplexUnitDisk{BigFloat,:open}()) == "ComplexUnitDisk{BigFloat}()  (open)"

        @test pseudolevel(ComplexUnitCircle(), 0.1) isa SublevelSet{Complex{Float64},:open}
        p = pseudolevel(ComplexUnitCircle(), 0.1)
        @test 0.8 ∉ p
        @test 0.95 ∈ p
        @test 1+0.1im ∈ p
        @test 1.1+0.2im ∉ p
    end

    @testset "spheres" begin
        # test UnitSphere constructor
        @test UnitSphere(2) isa VectorUnitSphere{Float64}
        @test UnitSphere(Val(2)) isa EuclideanUnitSphere{2,Float64}
        @test UnitSphere{Float64}() isa StaticUnitSphere{Float64}
        @test UnitSphere{Float64}(1) isa StaticUnitSphere{Float64}
        @test_throws AssertionError UnitSphere{Float64}(2)
        @test UnitSphere{SVector{2,Float64}}(Val(2)) isa EuclideanUnitSphere{2,Float64}
        @test_throws AssertionError UnitSphere{SVector{2,Float64}}(Val(3))
        @test UnitSphere{SVector{2,Float64}}(2) isa EuclideanUnitSphere{2,Float64}
        @test_throws AssertionError UnitSphere{SVector{2,Float64}}(3)
        @test UnitSphere{SVector{2,Float64}}() isa EuclideanUnitSphere{2,Float64}
        @test UnitSphere{Vector{Float64}}(2) isa VectorUnitSphere{Float64}
        @test_throws MethodError UnitSphere{Vector{Float64}}()

        @test DynamicUnitSphere{Float64}(1) isa DynamicUnitSphere{Float64}
        @test_throws AssertionError DynamicUnitSphere{Float64}(2)

        # and the generic sphere constructor
        @test Sphere() isa UnitSphere
        @test Sphere(1.0, 2.0) isa DomainSets.GenericSphere{Float64}
        @test Sphere{Float64}() isa UnitSphere{Float64}
        @test Sphere(1.0) == DomainSets.GenericSphere{SVector{3,Float64}}(1.0)
        @test Sphere{BigFloat}(1) == DomainSets.GenericSphere{BigFloat}(big(1))
        @test_throws MethodError Sphere{Vector{Float64}}(1.0)
        @test Sphere(1.0, [1,2,3]) isa DomainSets.GenericSphere{Vector{Float64},Float64}
        @test Sphere{Vector{Float64}}(1.0, [1,2,3]) isa DomainSets.GenericSphere{Vector{Float64},Float64}

        @test map_domain(AffineMap(2, [1,2,3]), UnitSphere()) isa DomainSets.GenericSphere
        @test map_domain(AffineMap(2, 3), UnitSphere{Float64}()) isa DomainSets.GenericSphere
        @test map_domain(LinearMap{SVector{3,Float64}}(2), UnitSphere()) isa DomainSets.GenericSphere
        @test map_domain(LinearMap(2), UnitSphere{Float64}()) isa DomainSets.GenericSphere
        @test map_domain(LinearMap(2), UnitSphere()) isa DomainSets.GenericSphere
        @test map_domain(Translation([1,2,3]), UnitSphere()) isa DomainSets.GenericSphere
        @test map_domain(Translation(1), UnitSphere{Float64}()) isa DomainSets.GenericSphere

        @test StaticUnitSphere() isa StaticUnitSphere{Float64}
        @test StaticUnitSphere(Val(2)) isa StaticUnitSphere{SVector{2,Float64}}

        @test GenericSphere() == GenericSphere(1.0)
        @test GenericSphere(2.0, 1:5) isa GenericSphere{Vector{Float64},Float64}
        @test GenericSphere(2, 1.0:5.0) isa GenericSphere{Vector{Float64},Float64}

        @test repr(UnitSphere()) == "UnitSphere()"
        @test repr(UnitSphere(Val(4))) == "UnitSphere(Val(4))"
        @test repr(UnitCircle()) == "UnitCircle()"
        @test repr(UnitCircle{BigFloat}()) == "UnitCircle{BigFloat}()"
        @test repr(VectorUnitSphere{Float64}(4)) == "UnitSphere(4)"
        @test repr(UnitSphere{Float64}()) == "UnitSphere{Float64}()"
        @test repr(Sphere(1.0,2.0)) == "Sphere(1.0, 2.0)"

        @test DomainSets.UnitCircleMap() == DomainSets.UnitCircleMap{Float64}()
        @test DomainSets.AngleMap() == DomainSets.AngleMap{Float64}()

        C = UnitCircle()
        @test SA[1.,0.] ∈ C
        @test SA[1.,1.] ∉ C
        @test approx_in(SA[1.,0.], C)
        @test !approx_in(SA[1.,1.], C)
        @test !isempty(C)
        @test isclosedset(C)
        @test !isopenset(C)
        @test parameterdomain(C) == UnitInterval()
        @test hasparameterization(C)
        p = parameterization(C)
        @test mapsize(p) == (2,)
        x = applymap(p, 1/2)
        @test jacobian(p, 0.4) ≈ SA[-2pi*sin(2pi*0.4), 2pi*cos(2pi*0.4)]
        @test diffvolume(p, 0.4) ≈ 2*pi
        @test diffvolume(p)(0.4) ≈ 2*pi
        @test approx_in(x, C)
        q = leftinverse(p)
        @test applymap(q, x) ≈ 1/2
        @test applymap(q, x) ≈ leftinverse(p, x)
        @test applymap(q, -x) ≈ 1
        @test rightinverse(q) == p
        @test rightinverse(q)(0) ≈ rightinverse(q, 0)
        @test jacobian(q) isa DomainSets.LazyJacobian
        @test jacobian(q, x) isa LinearAlgebra.Transpose{Float64,SVector{2,Float64}}

        @test boundingbox(C) == ProductDomain(ChebyshevInterval(), ChebyshevInterval())
        @test boundingbox(UnitSphere{Float64}()) == ChebyshevInterval()

        @test convert(LevelSet, UnitCircle()) isa LevelSet{SVector{2,Float64}}
        @test convert(LevelSet{SVector{2,BigFloat}}, UnitCircle()) isa LevelSet{SVector{2,BigFloat}}
        @test pseudolevel(UnitCircle(), 0.1) isa SublevelSet
        @test SA[1.05,0] ∈ pseudolevel(UnitCircle(), 0.1)
        @test SA[1.15,0] ∉ pseudolevel(UnitCircle(), 0.1)

        C2 = convert(Domain{SVector{2,BigFloat}}, C)
        @test eltype(C2) == SVector{2,BigFloat}

        @test convert(Domain{SVector{2,Float64}}, UnitSphere(2)) isa StaticUnitSphere
        @test convert(Domain{Vector{Float64}}, UnitSphere(Val(2))) isa DynamicUnitSphere

        C = 2UnitCircle() .+ SA[1.,1.]
        @test C isa DomainSets.GenericSphere
        @test approx_in(SA[3.,1.], C)

        @test canonicaldomain(C) == UnitCircle()
        @test matrix(mapfrom_canonical(C)) == [2 0; 0 2]
        @test vector(mapfrom_canonical(C)) == [1; 1]
        @test parameterdomain(C) == UnitInterval()
        @test mapfrom_parameterdomain(C) isa ComposedMap
        @test mapfrom_parameterdomain(C)(0.5) ≈ [-1; 1]
        @test boundingbox(C) == (-1.0..3.0)^2

        C = UnitCircle() .+ SA[1,1]
        @test C isa DomainSets.GenericSphere
        @test approx_in(SA[2,1], C)

        S = UnitSphere()
        @test SA[1.,0.,0.] ∈ S
        @test SA[1.,0.,1.] ∉ S
        @test approx_in(SA[cos(1.),sin(1.),0.], S)
        @test !isempty(S)
        S2 = convert(Domain{SVector{3,BigFloat}}, S)
        @test eltype(S2) == SVector{3,BigFloat}

        @test Basis3Vector() in S

        @test issubset(UnitSphere(), UnitBall())

        S = 2 * UnitSphere() .+ SA[1.,1.,1.]
        @test S isa DomainSets.GenericSphere
        @test approx_in(SA[1. + 2*cos(1.),1. + 2*sin(1.),1.], S)
        @test !approx_in(SA[4.,1.,5.], S)

        D = UnitCircle()
        @test convert(Domain{SVector{2,BigFloat}}, D) ≡ UnitCircle{BigFloat}()
        @test SVector(1,0) in D
        @test SVector(nextfloat(1.0),0) ∉ D

        D = UnitSphere()
        @test convert(Domain{SVector{3,BigFloat}}, D) ≡ UnitSphere{SVector{3,BigFloat}}()
        @test SVector(1,0,0) in D
        @test SVector(nextfloat(1.0),0,0) ∉ D

        D4 = UnitSphere(4)
        @test D4 isa DynamicUnitSphere
        @test D4 isa VectorUnitSphere
        @test dimension(D4) == 4
        cheb = ChebyshevInterval()
        @test boundingbox(D4) == ProductDomain([cheb, cheb, cheb, cheb])
    end

    @testset "derived types" begin
        # Create an ellipse, the curve
        E = ellipse(2.0, 4.0)
        @test SA[2.0,0.0] ∈ E
        @test SA[0.0,4.0] ∈ E
        @test SA[2.0+1e-10,0.0] ∉ E
        @test SA[0.0,0.0] ∉ E
        E = ellipse(1, 2.0)
        @test eltype(E) == SVector{2,Float64}

        # Create an ellipse, the domain with volume
        E2 = ellipse_shape(2.0, 4.0)
        @test SA[2.0,0.0] ∈ E2
        @test SA[0.0,4.0] ∈ E2
        @test SA[2.0+1e-10,0.0] ∉ E2
        @test SA[0.0,0.0] ∈ E2
        @test SA[1.0,1.0] ∈ E2

        E2 = ellipse_shape(1, 2.0)
        @test eltype(E) == SVector{2,Float64}

        C = cylinder()
        @test eltype(C) == SVector{3,Float64}
        C2 = cylinder(1.0, 2)
        @test SA[0.5,0.2,1.5] ∈ C2
    end

    @testset "mapped_domain" begin
        @test MappedDomain(cos, 0..1.0) isa MappedDomain{Float64}
        @test MappedDomain{Float64}(cos, 0..1.0) isa MappedDomain{Float64}
        @test cos.(0..1.0) isa MappedDomain
        @test isempty(MappedDomain(LinearMap(2.0), EmptySpace()))

        # Test chaining of maps
        D = UnitCircle()
        D1 = MappedDomain(inverse(LinearMap(2)), D)
        @test typeof(D1) <: MappedDomain
        @test typeof(superdomain(D1)) <: UnitSphere
        @test isclosedset(D1)
        @test !isopenset(D1)
        @test convert(Domain{SVector{2,BigFloat}}, D1) isa MappedDomain{SVector{2,BigFloat}}
        D2 = 2 * D1
        @test typeof(superdomain(D2)) <: UnitSphere

        D = UnitInterval()^2
        show(io, textmime, rotate(D,1.))
        @test String(take!(io))[1:17] == "A .* UnitSquare()"

        D = rotate(UnitInterval()^2, π)
        @test SA[-0.9, -0.9] ∈ D
        @test SA[-1.1, -1.1] ∉ D
        x = point_in_domain(D)
        @test forward_map(D)(x) ≈ forward_map(D, x)
        @test DomainSets.toexternalpoint(D, x) ≈ forward_map(D, x)

        D = rotate(UnitInterval()^2, π, SA[-.5,-.5])
        @test SA[-1.5, -1.5] ∈ D
        @test SA[-0.5, -0.5] ∉ D

        D = rotate(UnitInterval()^3 .+ SA[-.5,-.5,-.5], pi, pi, pi)
        @test SA[0.4, 0.4, 0.4] ∈ D
        @test SA[0.6, 0.6, 0.6] ∉ D

        D = rotate((-1.5.. 2.2) × (0.5 .. 0.7) × (-3.0 .. -1.0), π, π, π, SA[.35, .65, -2.])
        @test SA[0.9, 0.6, -2.5] ∈ D
        @test SA[0.0, 0.6, 0.0] ∉ D

        B = mapped_domain(inverse(LinearMap(2.0)), VectorUnitBall(10))
        @test dimension(B) == 10
        @test superdomain(B) ∘ inverse_map(B) == B
        @test isopenset(interior(B))
        @test B == closure(interior(B))
        @test DomainSets.superdomain(boundary(B)) isa UnitSphere
        @test canonicaldomain(B) == VectorUnitBall(10)
        @test mapfrom_canonical(B) == forward_map(B)
        @test mapto_canonical(B) == inverse_map(B)
        @test parameterdomain(B) == canonicaldomain(B)
        @test mapfrom_parameterdomain(B) == mapfrom_canonical(B)
        @test mapto_parameterdomain(B) == mapto_canonical(B)

        # Test parametric domain
        using DomainSets: ParametricDomain
        m = AffineMap(ones(2), [4; 5])
        pd = ParametricDomain(m, UnitInterval())
        @test pd isa Domain{Vector{Float64}}
        @test forward_map(pd) == m
        @test forward_map(pd, 0.4) ≈ m(0.4)
        @test mapfrom_canonical(pd) == m
        @test canonicaldomain(pd) == 0..1
        @test boundary(pd) isa UnionDomain
        @test interior(pd) isa ParametricDomain
        @test closure(pd) isa ParametricDomain
    end

    @testset "simplex" begin
        @test UnitSimplex(2) isa VectorUnitSimplex{Float64}
        @test UnitSimplex(Val(2)) isa EuclideanUnitSimplex{2,Float64}
        @test UnitSimplex{Float64}() isa StaticUnitSimplex{Float64}
        @test UnitSimplex{Float64}(1) isa StaticUnitSimplex{Float64}
        @test_throws AssertionError UnitSimplex{Float64}(2)
        @test UnitSimplex{SVector{2,Float64}}(Val(2)) isa EuclideanUnitSimplex{2,Float64}
        @test_throws AssertionError UnitSimplex{SVector{2,Float64}}(Val(3))
        @test UnitSimplex{SVector{2,Float64}}(2) isa EuclideanUnitSimplex{2,Float64}
        @test_throws AssertionError UnitSimplex{SVector{2,Float64}}(3)
        @test UnitSimplex{SVector{2,Float64}}() isa EuclideanUnitSimplex{2,Float64}
        @test UnitSimplex{Vector{Float64}}(2) isa VectorUnitSimplex{Float64}
        @test_throws MethodError UnitSimplex{Vector{Float64}}()

        @test UnitSimplex{Float64,:open}() isa StaticUnitSimplex{Float64,:open}
        @test UnitSimplex{Float64,:closed}(1) isa StaticUnitSimplex{Float64,:closed}
        @test_throws AssertionError UnitSimplex{Float64,:closed}(2)
        @test UnitSimplex{SVector{2,Float64},:open}(Val(2)) isa EuclideanUnitSimplex{2,Float64,:open}
        @test_throws AssertionError UnitSimplex{SVector{2,Float64},:closed}(Val(3))
        @test UnitSimplex{SVector{2,Float64},:open}(2) isa EuclideanUnitSimplex{2,Float64,:open}
        @test_throws AssertionError UnitSimplex{SVector{2,Float64},:closed}(3)
        @test UnitSimplex{SVector{2,Float64},:open}() isa EuclideanUnitSimplex{2,Float64,:open}
        @test UnitSimplex{Vector{Float64},:closed}(2) isa VectorUnitSimplex{Float64,:closed}
        @test_throws MethodError UnitSimplex{Vector{Float64},:open}()

        @test StaticUnitSimplex(Val(3)) isa StaticUnitSimplex{SVector{3,Float64}}
        @test DynamicUnitSimplex{Float64}(1) isa DynamicUnitSimplex{Float64}
        @test_throws AssertionError DynamicUnitSimplex{Float64}(2)

        @test repr(UnitSimplex(Val(2))) == "UnitSimplex(Val(2))"
        @test repr(UnitSimplex(3)) == "UnitSimplex(3)"

        d = UnitSimplex(Val(2))
        # We test a point in the interior, a point on each of the boundaries and
        # all corners.
        @test SA[0.2,0.2] ∈ d
        @test SA[0.0,0.2] ∈ d
        @test SA[0.2,0.0] ∈ d
        @test SA[0.5,0.5] ∈ d
        @test SA[0.0,0.0] ∈ d
        @test SA[1.0,0.0] ∈ d
        @test SA[0.0,1.0] ∈ d
        # And then some points outside
        @test SA[0.6,0.5] ∉ d
        @test SA[0.5,0.6] ∉ d
        @test SA[-0.2,0.2] ∉ d
        @test SA[0.2,-0.2] ∉ d
        @test boundingbox(d) == UnitCube{SVector{2,Float64}}()

        @test approx_in(SA[-0.1,-0.1], d, 0.1)
        @test !approx_in(SA[-0.1,-0.1], d, 0.09)

        @test corners(d) == [ SA[0.0,0.0], SA[1.0,0.0], SA[0.0,1.0]]

        @test convert(Domain{SVector{2,BigFloat}}, d) == EuclideanUnitSimplex{2,BigFloat}()

        @test isclosedset(d)
        @test !isopenset(d)
        @test isopenset(interior(d))
        @test closure(d) == d
        @test point_in_domain(d) ∈ d

        # open/closed
        d2 = EuclideanUnitSimplex{2,Float64,:open}()
        @test !isclosedset(d2)
        @test isopenset(d2)
        @test SA[0.3,0.1] ∈ d2
        @test SA[0.0,0.1] ∉ d2
        @test SA[0.3,0.0] ∉ d2
        @test approx_in(SA[-0.01,0.0], d2, 0.1)
        @test !approx_in(SA[-0.01,0.0], d2, 0.001)

        d3 = EuclideanUnitSimplex{3,BigFloat}()
        @test point_in_domain(d3) ∈ d3
        x0 = big(0.0)
        x1 = big(1.0)
        x2 = big(0.3)
        @test SA[x0,x0,x0] ∈ d3
        @test SA[x1,x0,x0] ∈ d3
        @test SA[x0,x1,x0] ∈ d3
        @test SA[x0,x0,x1] ∈ d3
        @test SA[x2,x0,x0] ∈ d3
        @test SA[x0,x2,x0] ∈ d3
        @test SA[x0,x0,x2] ∈ d3
        @test SA[x2,x2,x2] ∈ d3
        @test SA[-x2,x2,x2] ∉ d3
        @test SA[x2,-x2,x2] ∉ d3
        @test SA[x2,x2,-x2] ∉ d3
        @test SA[x1,x1,x1] ∉ d3

        D = VectorUnitSimplex(2)
        @test isopenset(interior(D))
        @test closure(D) == D
        @test SA[0.2,0.2] ∈ D
        @test SA[0.0,0.2] ∈ D
        @test SA[0.2,0.0] ∈ D
        @test SA[0.5,0.5] ∈ D
        @test SA[0.0,0.0] ∈ D
        @test SA[1.0,0.0] ∈ D
        @test SA[0.0,1.0] ∈ D
        # And then some points outside
        @test SA[0.6,0.5] ∉ D
        @test SA[0.5,0.6] ∉ D
        @test SA[-0.2,0.2] ∉ D
        @test SA[0.2,-0.2] ∉ D
        @test convert(Domain{Vector{BigFloat}}, D) == VectorUnitSimplex{BigFloat}(2)
        @test corners(D) == [ [0.0,0.0], [1.0,0.0], [0.0,1.0]]
        @test boundingbox(D) == UnitCube(4)
    end

    @testset "level sets" begin
        d1 = LevelSet(cos, 1.0)
        @test d1 isa LevelSet{Float64}
        @test convert(Domain{ComplexF64}, d1) isa LevelSet{ComplexF64}
        show(io,d1)
        @test String(take!(io)) == "level set f(x) = 1.0 with f = cos"
        @test 0.0 ∈ d1
        @test 0im ∈ d1
        @test 0.1 ∉ d1
        @test 0.1+1im ∉ d1

        # prod yields the function (x,y) -> x*y
        d2 = ZeroSet{SVector{2,Float64}}(prod)
        @test d2 isa ZeroSet{SVector{2,Float64}}
        @test SA[0.1,0.3] ∉ d2
        @test SA[0.0,0.3] ∈ d2
        @test SA[0.1,0.0] ∈ d2
        @test ZeroSet(cos) isa ZeroSet{Float64}
        @test convert(Domain{BigFloat}, ZeroSet(cos)) isa ZeroSet{BigFloat}
        @test convert(LevelSet, ZeroSet{BigFloat}(cos)) isa LevelSet{BigFloat}
        @test convert(LevelSet{BigFloat}, ZeroSet{Float64}(cos)) isa LevelSet{BigFloat}

        d3 = SublevelSet(cos, 0.5)
        d3_open = SublevelSet{Float64,:open}(cos,0.5)
        @test d3 isa SublevelSet{Float64,:closed}
        @test interior(d3) == d3_open
        @test closure(d3_open) == d3
        @test closure(d3) == d3
        @test interior(d3_open) == d3_open
        @test boundary(d3) == LevelSet(cos, 0.5)
        @test 3.0 ∈ d3
        @test 0.0 ∉ d3
        @test 0.0 ∉ d3_open
        show(io, d3)
        @test String(take!(io)) == "sublevel set f(x) <= 0.5 with f = cos"
        show(io, d3_open)
        @test String(take!(io)) == "sublevel set f(x) < 0.5 with f = cos"
        @test convert(Domain{BigFloat}, d3) isa SublevelSet{BigFloat,:closed}
        @test convert(Domain{BigFloat}, d3_open) isa SublevelSet{BigFloat,:open}

        d4 = SubzeroSet{SVector{2,Float64}}(prod)
        d4_open = SubzeroSet{SVector{2,Float64},:open}(prod)
        @test d4 isa SubzeroSet{SVector{2,Float64},:closed}
        @test interior(d4) == d4_open
        @test closure(d4_open) == d4
        @test closure(d4) == d4
        @test interior(d4_open) == d4_open
        @test boundary(d4) == ZeroSet{SVector{2,Float64}}(prod)
        @test SA[0.1,0.3] ∉ d4
        @test SA[-0.1,0.3] ∈ d4
        @test SA[-0.1,-0.3] ∉ d4
        @test SA[-0.1,0.3] ∈ d4_open
        convert(Domain{SVector{2,BigFloat}}, d4) isa SubzeroSet{SVector{2,BigFloat},:closed}
        convert(Domain{SVector{2,BigFloat}}, d4_open) isa SubzeroSet{SVector{2,BigFloat},:open}
        @test SubzeroSet(cos) == SubzeroSet{Float64}(cos)

        d5 = SuperlevelSet(cos, 0.5)
        d5_open = SuperlevelSet{Float64,:open}(cos, 0.5)
        @test d5 isa SuperlevelSet{Float64,:closed}
        @test interior(d5) == d5_open
        @test closure(d5_open) == d5
        @test closure(d5) == d5
        @test interior(d5_open) == d5_open
        @test boundary(d5) == LevelSet(cos, 0.5)
        @test 3.0 ∉ d5
        @test 0.0 ∈ d5
        @test 0.0 ∈ d5
        @test 0.0 ∈ d5_open
        show(io, d5)
        @test String(take!(io)) == "superlevel set f(x) >= 0.5 with f = cos"
        show(io, d5_open)
        @test String(take!(io)) == "superlevel set f(x) > 0.5 with f = cos"
        @test convert(Domain{BigFloat}, d5) isa SuperlevelSet{BigFloat}
        @test convert(Domain{BigFloat}, d5_open) isa SuperlevelSet{BigFloat,:open}

        d6 = SuperzeroSet{SVector{2,Float64}}(prod)
        d6_open = SuperzeroSet{SVector{2,Float64},:open}(prod)
        @test d6 isa SuperzeroSet{SVector{2,Float64},:closed}
        @test interior(d6) == d6_open
        @test closure(d6_open) == d6
        @test closure(d6) == d6
        @test interior(d6_open) == d6_open
        @test boundary(d6) == ZeroSet{SVector{2,Float64}}(prod)
        @test SA[0.1,0.3] ∈ d6
        @test SA[-0.1,0.3] ∉ d6
        @test SA[-0.1,-0.3] ∈ d6
        @test SuperzeroSet(cos) isa SuperzeroSet{Float64}
        @test convert(Domain{SVector{2,BigFloat}}, d6) isa SuperzeroSet{SVector{2,BigFloat},:closed}
        @test convert(Domain{SVector{2,BigFloat}}, d6_open) isa SuperzeroSet{SVector{2,BigFloat},:open}
    end

    @testset "indicator functions" begin
        ispositive(x) = x >= 0
        d = IndicatorFunction(ispositive)
        @test d isa IndicatorFunction{Float64}
        @test DomainSets.indicatorfunction(d) == ispositive
        show(io,d)
        @test String(take!(io)) == "indicator domain defined by function f = ispositive"
        @test 0 ∈ d
        @test big(0) ∈ d
        @test -1 ∉ d
        @test d ∩ ChebyshevInterval() isa DomainSets.BoundedIndicatorFunction
        @test ChebyshevInterval() ∩ d isa DomainSets.BoundedIndicatorFunction

        @test convert(IndicatorFunction, 0..1) isa IndicatorFunction
        @test convert(IndicatorFunction, d) == d
        @test convert(Domain{BigFloat}, d) isa IndicatorFunction{BigFloat}
        @test 0.5 ∈ convert(IndicatorFunction, 0..1)

        d2 = Domain(x>0 for x in -1..1)
        @test -0.5 ∉ d2
        @test 0.5 ∈ d2
        show(io, d2)
        @test String(take!(io)) == "indicator function bounded by: -1..1"

        d3 = Domain(x*y>0 for (x,y) in UnitDisk())
        @test [0.4,0.2] ∈ d3
        @test [0.4,-0.2] ∉ d3

        d4 = Domain( x+y+z > 0 for (x,y) in UnitDisk(), z in 0..1)
        @test d4 isa DomainSets.BoundedIndicatorFunction{F,<:TupleProductDomain} where F
        @test DomainSets.indicatorfunction(d4) isa Function
        @test ( [0.5,0.2], 0.5) ∈ d4
        @test ( [0.5,0.2], 1.5) ∉ d4
        @test ( [-0.5,-0.2], 0.1) ∉ d4
        @test boundingbox(d4) == boundingbox(DomainSets.boundingdomain(d4))
    end
end

@testset "cartesian product" begin
    @test productdomain() == ()
    @test productdomain(2) == 2

    @testset "VcatDomain" begin
        d1 = VcatDomain(-1.0..1.0, -1.0..1.0)
        @test d1 isa DomainSets.VcatDomain
        @test d1.domains isa Tuple
        @test eltype(d1) == SVector{2,typeof(1.0)}
        @test SA[0.5,0.5] ∈ d1
        @test SA[-1.1,0.3] ∉ d1
        @test @inferred(VcatDomain(-1.0..1, -1.0..1)) === d1
        # Test promotion
        @test convert(Domain{SVector{2,BigFloat}}, d1) isa ProductDomain{SVector{2,BigFloat}}
        d1w = convert(Domain{SVector{2,BigFloat}}, d1)
        @test eltype(d1w) == SVector{2,BigFloat}
        @test eltype(component(d1w, 1)) == BigFloat
        @test eltype(component(d1w, 2)) == BigFloat

        @test VcatDomain( (-1..1, -2..2)) isa VcatDomain{2,Int,(1,1),Tuple{ClosedInterval{Int64}, ClosedInterval{Int64}}}

        show(io,d1)
        @test String(take!(io)) == "(-1.0..1.0) × (-1.0..1.0)"

        bnd = boundary(d1)
        @test bnd isa EuclideanDomain
        @test bnd isa UnionDomain
        @test [-1.0, 0.5] ∈ bnd
        @test [1.0, 0.5] ∈ bnd
        @test [0.5, -1.0] ∈ bnd
        @test [0.5, 1.0] ∈ bnd
        @test [0.5, 0.2] ∉ bnd

        # A single domain
        @test VcatDomain(UnitCircle()) isa VcatDomain{2}

        # Test vectors of wrong length
        @test_logs (:warn, "`in`: incompatible combination of vector with length 3 and domain '(-1.0..1.0) × (-1.0..1.0)' with dimension 2. Returning false.") SA[0.0,0.0,0.0] ∉ d1
        @test_logs (:warn, "`in`: incompatible combination of vector with length 1 and domain '(-1.0..1.0) × (-1.0..1.0)' with dimension 2. Returning false.") SA[0.0] ∉ d1
        @test_logs (:warn, "`in`: incompatible combination of vector with length 3 and domain '(-1.0..1.0) × (-1.0..1.0)' with dimension 2. Returning false.") [0.0,0.0,0.0] ∉ d1
        @test_logs (:warn, "`in`: incompatible combination of vector with length 1 and domain '(-1.0..1.0) × (-1.0..1.0)' with dimension 2. Returning false.") [0.0] ∉ d1

        d3 = VcatDomain(-1.0 .. 1.0, -1.5 .. 2.5)
        @test SA[0.5,0.5] ∈ d3
        @test SA[-1.1,0.3] ∉ d3

        d3 = VcatDomain(1.05 * UnitDisk(), -1.0 .. 1.0)
        @inferred(cross(1.05 * UnitDisk(), -1.0 .. 1.0)) === d3
        @test d3 isa VcatDomain
        @test eltype(d3) == SVector{3,Float64}
        @test SA[0.5,0.5,0.8] ∈ d3
        @test SA[-1.1,0.3,0.1] ∉ d3
        @test point_in_domain(d3) ∈ d3
    end
    @testset "mixed intervals" begin
        d = (0..1) × (0.0..1)
        @test SA[0.1,0.2] ∈ d
        @test SA[0.1,1.2] ∉ d
        @test SA[1.1,1.3] ∉ d
        @test d isa EuclideanDomain{2}
        # Make sure promotion of domains happened
        @test eltype(component(d,1)) == Float64
        @test eltype(component(d,2)) == Float64
        @test point_in_domain(d) ∈ d
    end
    @testset "vector domains" begin
        d1 = VectorProductDomain([0..1.0, 0..2.0])
        @test d1 isa VectorDomain{Float64}
        @test d1.domains isa Vector
        @test dimension(d1) == 2
        @test [0.1,0.2] ∈ d1
        @test SA[0.1,0.2] ∈ d1
        @test point_in_domain(d1) ∈ d1
        @test convert(Domain{Vector{BigFloat}}, d1) == d1
        d1big = convert(Domain{Vector{BigFloat}}, d1)
        @test eltype(d1big) == Vector{BigFloat}

        # Test an integer type as well
        d2 = VectorProductDomain([0..1, 0..3])
        @test dimension(d2) == 2
        @test [0.1,0.2] ∈ d2
        @test point_in_domain(d2) ∈ d2

        # other constructor calls
        @test VectorProductDomain(0..1, 1..2) isa VectorProductDomain{Vector{Int}}
        @test VectorProductDomain((0..1, 1..2)) isa VectorProductDomain{Vector{Int}}
        @test VectorProductDomain{Vector{Int}}(0..1, 1..2) isa VectorProductDomain{Vector{Int}}
        @test VectorProductDomain{Vector{Int}}((0..1, 1..2)) isa VectorProductDomain{Vector{Int}}

        bnd = boundary(d1)
        @test bnd isa VectorDomain
        @test bnd isa UnionDomain
        @test dimension(bnd) == 2
        @test [0.0, 0.5] ∈ bnd
        @test [1.0, 0.5] ∈ bnd
        @test [0.2, 0.0] ∈ bnd
        @test [0.2, 2.0] ∈ bnd
        @test [0.2, 0.5] ∉ bnd

        d3 = VectorProductDomain([0..i for i in 1:10])
        @test d3 isa VectorProductDomain
        @test eltype(d3) == Vector{Int}
        @test rand(10) ∈ d3
        @test 2 .+ rand(10) ∉ d3
        d4 = VectorProductDomain{Vector{Float64}}([0..i for i in 1:10])
        @test d4 isa VectorProductDomain
        @test eltype(d4) == Vector{Float64}
        @test rand(10) ∈ d4
        @test 2 .+ rand(10) ∉ d4

        @test VectorProductDomain{SVector{2,Float64}}(SVector(0..1, 0..2)).domains[1] isa Domain{Float64}
    end
    @testset "Tuple product domain" begin
        # Use the constructor ProductDomain{T} directly
        d1 = TupleProductDomain{Tuple{Float64,Float64}}(0..0.5, 0..0.7)
        @test d1 isa TupleProductDomain{Tuple{Float64,Float64}}
        @test d1.domains isa Tuple
        @test eltype(d1) == Tuple{Float64,Float64}
        @test dimension(d1) == 2
        @test (0.2,0.6) ∈ d1
        @test (0.2,0.8) ∉ d1
        @test (true,0.6) ∉ d1
        if VERSION < v"1.6-"
            @test_logs (:warn, "`in`: incompatible combination of point: SArray{Tuple{2},Float64,1,2} and domain eltype: Tuple{Float64,Float64}. Returning false.") SA[0.2,0.6] ∉ d1
            @test_logs (:warn, "`in`: incompatible combination of point: Array{Float64,1} and domain eltype: Tuple{Float64,Float64}. Returning false.") [0.2,0.6] ∉ d1
        else
            @test_logs (:warn, "`in`: incompatible combination of point: SVector{2, Float64} and domain eltype: Tuple{Float64, Float64}. Returning false.") SA[0.2,0.6] ∉ d1
            @test_logs (:warn, "`in`: incompatible combination of point: Vector{Float64} and domain eltype: Tuple{Float64, Float64}. Returning false.") [0.2,0.6] ∉ d1
        end
        @test convert(Domain{Tuple{BigFloat,BigFloat}}, d1) == d1
        d1big = convert(Domain{Tuple{BigFloat,BigFloat}}, d1)
        @test eltype(d1big) == Tuple{BigFloat,BigFloat}
        @test eltype(component(d1big,1)) == BigFloat
        @test eltype(component(d1big,2)) == BigFloat

        d2 = TupleProductDomain(['a','b'], 0..1)
        @test d2.domains isa Tuple
        @test dimension(d2) == 2
        @test eltype(d2) == Tuple{Char,Int}
        @test ('a',0.4) ∈ d2
        @test ('b',1.5) ∉ d2
        @test ('c',0.5) ∉ d2

        # other constructor calls
        @test TupleProductDomain(0..1, 1..2.0) isa TupleProductDomain{Tuple{Int,Float64}}
        @test TupleProductDomain([0..1, 1..2.0]) isa TupleProductDomain{Tuple{Float64,Float64}}
        @test TupleProductDomain{Tuple{Int,Float64}}(0..1, 1..2.0) isa TupleProductDomain{Tuple{Int,Float64}}
        @test TupleProductDomain{Tuple{Float64,Float64}}(0..1, 1..2.0) isa TupleProductDomain{Tuple{Float64,Float64}}
        @test TupleProductDomain{Tuple{Float64,Float64}}([0..1, 1..2.0]) isa TupleProductDomain{Tuple{Float64,Float64}}

        bnd = boundary(d1)
        @test eltype(bnd) == Tuple{Float64,Float64}
        @test bnd isa UnionDomain
        @test dimension(bnd) == 2
        @test (0.0, 0.5) ∈ bnd
        @test (0.5, 0.5) ∈ bnd
        @test (0.3, 0.2) ∉ bnd
        @test (0.3, 0.0) ∈ bnd
        @test (0.3, 0.7) ∈ bnd
    end
    @testset "cube" begin
        @test volume(UnitCube()) == 1
        @test EuclideanUnitCube{2}() == EuclideanUnitCube{2,Float64}()
        @test UnitSquare() == UnitSquare{Float64}()
        @test UnitCube(Val(2)) isa EuclideanUnitCube{2,Float64}
        @test UnitCube{SVector{2,BigFloat}}(Val(2)) isa EuclideanUnitCube{2,BigFloat}

        @test Set(corners((0..1)^2)) == Set([ [1,0], [0,1], [1,1], [0,0]])
        @test Set(corners(UnitCube())) == Set([ [0,0,0], [1,0,0], [0,1,0], [0,0,1], [1,1,0], [1,0,1], [0,1,1], [1,1,1]])
        @test boundary(boundary(UnitSquare())) == UnionDomain(map(t->Point(SVector{2}(t)), corners(UnitSquare())))
        @test boundary(boundary(boundary(UnitCube()))) == UnionDomain(map(t->Point(SVector{3}(t)), corners(UnitCube())))
        @test UnitCube() == UnitCube(3)
        @test boundary(UnitCube()) == boundary(UnitCube(3))
        @test boundary(boundary(UnitCube())) == boundary(boundary(UnitCube(3)))
        @test boundary(boundary(boundary(UnitCube()))) == boundary(boundary(boundary(UnitCube(3))))
        @test UnitSquare() == UnitCube(2)
        @test boundary(UnitSquare()) == boundary(UnitCube(2))
        @test boundary(boundary(UnitSquare())) == boundary(boundary(UnitCube(2)))

        @test StaticUnitCube() == UnitCube()
        @test StaticUnitCube(Val(3)) == UnitCube()
        @test convert(Domain{SVector{2,Float64}}, UnitCube()) == StaticUnitCube(Val(2))
        @test convert(Domain{Vector{Float64}}, UnitCube()) === VectorUnitCube(3)

        @test DynamicUnitCube{SVector{3,Float64}}(3) isa DynamicUnitCube
        @test_throws AssertionError DynamicUnitCube{SVector{3,Float64}}(2)

        d1 = VectorUnitCube{Float64}(4)
        @test VectorUnitCube(4) == d1
        @test UnitCube(4) == d1
        @test UnitCube{Vector{Float64}}(4) == d1
        @test dimension(d1) == 4
        @test component(d1, 1) == 0..1
        @test SA[0.9,0.9,0.4,0.2] ∈ d1
        @test [1.2,0.3,0.4,0.6] ∉ d1
        @test convert(Domain{SVector{4,Float64}}, d1) isa StaticUnitCube
        @test convert(Domain{Vector{BigFloat}}, d1) === VectorUnitCube{BigFloat}(4)
        @test DomainSets.VectorUnitSquare() === UnitCube(2)

        @test ProductDomain([UnitInterval(),UnitInterval()]) isa VectorUnitCube{Float64}
        @test ProductDomain{Vector{Float64}}([UnitInterval(),UnitInterval()]) isa VectorUnitCube{Float64}
        @test ProductDomain{Vector{BigFloat}}([UnitInterval(),UnitInterval()]) isa VectorUnitCube{BigFloat}
        @test ProductDomain{SVector{2,BigFloat}}(UnitInterval(),UnitInterval()) isa EuclideanUnitCube{2,BigFloat}
        @test ProductDomain{SVector{2,BigFloat}}((UnitInterval(),UnitInterval())) isa EuclideanUnitCube{2,BigFloat}
        @test ProductDomain{SVector{2,BigFloat}}(SVector(UnitInterval(),UnitInterval())) isa EuclideanUnitCube{2,BigFloat}

        @test repr(UnitCube()) == "UnitCube()"
        @test repr(UnitCube(Val(4))) == "UnitCube(Val(4))"
        @test repr(UnitSquare()) == "UnitSquare()"
        @test repr(UnitSquare{BigFloat}()) == "UnitSquare{BigFloat}()"
        @test repr(UnitCube(4)) == "UnitCube(4)"

        D = UnitInterval()^2
        @test D isa EuclideanUnitCube{2,Float64}
        @test SA[0.9, 0.9] ∈ D
        @test SA[1.1, 1.1] ∉ D
        @test !isempty(D)
        @test isclosedset(D)
        @test !isopenset(D)
        @test convert(DomainSets.HyperRectangle, VcatDomain(UnitInterval(), UnitInterval())) isa StaticUnitCube

        @test approx_in(SA[-0.1,-0.1], D, 0.1)
        @test !approx_in(SA[-0.1,-0.1], D, 0.09)

        @test intersect(ChebyshevInterval()^3, UnitInterval()^3) isa EuclideanUnitCube{3,Float64}

        #Cube
        D = (-1.5 .. 2.2) × (0.5 .. 0.7) × (-3.0 .. -1.0)
        @test SA[0.9, 0.6, -2.5] ∈ D
        @test SA[0.0, 0.6, 0.0] ∉ D

        @test issubset( (0..1)^3, (-1..2)^3 )
    end
    @testset "Rectangle" begin
        d1 = (-1.0..1.0) × (-1.0..1.0)

        d4 = d1 × (-1.0..1.0)
        @test d4 isa Rectangle
        @test SA[0.5,0.5,0.8] ∈ d4
        @test SA[-1.1,0.3,0.1] ∉ d4
        @test point_in_domain(d4) ∈ d4

        @test d1[Component(1)] == -1..1
        @test d1[Component(2)] == -1..1
        @test_throws BoundsError d1[Component(3)]

        d5 = (-1.0..1.)×d1
        @test d5 isa Rectangle
        @test SA[0.,0.5,0.5] ∈ d5
        @test SA[0.,-1.1,0.3] ∉ d5
        @test point_in_domain(d5) ∈ d5

        d6 = d1 × d1
        @test d6 isa Rectangle
        @test SA[0.,0.,0.5,0.5] ∈ d6
        @test SA[0.,0.,-1.1,0.3] ∉ d6
        @test point_in_domain(d6) ∈ d6

        @test Rectangle( SA[1,2], SA[2.0,3.0]) isa Rectangle{SVector{2,Float64}}
        @test Rectangle([0..1, 2..3]) isa Rectangle{Vector{Int}}
        @test Rectangle((0..1, 2..3)) isa Rectangle{SVector{2,Int}}
        @test Rectangle{SVector{2,Float64}}((0..1, 2..3)) isa Rectangle{SVector{2,Float64}}

        @test_throws ErrorException Rectangle(UnitCircle(), UnitDisk())
        @test_throws ErrorException Rectangle(OpenInterval(1,2), 3..4)
        @test_throws ErrorException Rectangle{SVector{2,Float64}}(UnitCircle(), UnitDisk())

        bnd = boundary(Rectangle([1,2],[3,4]))
        @test [1,3] ∈ bnd
        @test [1,2.5] ∈ bnd
        @test [1.5,4] ∈ bnd
        @test [1.5,3.5] ∉ bnd
    end
    @testset "fixed product domains" begin
        d1 = ProductDomain(ChebyshevInterval(), ChebyshevInterval())
        @test d1 isa DomainSets.FixedIntervalProduct
        @test d1 isa DomainSets.ChebyshevProductDomain
        @test volume(d1) == 4
        @test component(d1,1) isa ChebyshevInterval{Float64}
        @test components(d1) isa NTuple{2,ChebyshevInterval{Float64}}
        @test convert(Domain{SVector{2,BigFloat}}, d1) isa DomainSets.ChebyshevProductDomain{2,BigFloat}
        @test d1 == DomainSets.ChebyshevProductDomain(Val(2))
        @test d1 == DomainSets.ChebyshevProductDomain{2}()
        @test ProductDomain(components(d1)) == d1
        @test ProductDomain{SVector{2,BigFloat}}(components(d1)) isa DomainSets.ChebyshevProductDomain{2,BigFloat}
        @test ProductDomain{SVector{2,BigFloat}}(components(d1)...) isa DomainSets.ChebyshevProductDomain{2,BigFloat}
    end
    @testset "ProductDomain" begin
        d1 = 0..1.0
        d2 = UnitBall{Float64}()
        d3 = UnitCircle()
        @test ProductDomain{SVector{2,Float64}}(d1, d2) isa VcatDomain
        @test ProductDomain{Tuple{Float64,Float64}}(d1, d2) isa TupleProductDomain
        @test ProductDomain{Vector{Float64}}([d1; d2]) isa VectorProductDomain

        @test ProductDomain((d1,d3)) isa VcatDomain
        @test ProductDomain((d1,d1)) isa Rectangle

        @test volume(ProductDomain(0..1,1..3)) == 2

        @test ProductDomain(SVector(0..1, 0..2)) isa Rectangle{SVector{2,Int}}
        @test ProductDomain(1.05 * UnitDisk(), -1.0 .. 1.0) isa VcatDomain{3,Float64}
        @test ProductDomain(['a','b'], 0..1) isa TupleProductDomain

        @test ProductDomain{Tuple{Float64,Float64}}(0..0.5, 0..0.7) isa Rectangle{Tuple{Float64,Float64}}
        # Some conversions
        @test convert(Domain{Vector{Float64}}, productdomain(d1,d2)) isa VectorProductDomain{Vector{Float64}}
        @test convert(Domain{SVector{2,Float64}}, ProductDomain([d1,d2])) isa VcatDomain{2,Float64}
        @test convert(Domain{Vector{Float64}}, TupleProductDomain(d1,d2)) isa VectorDomain{Float64}

        # intersection of product domains
        @test ProductDomain([0..1.0, 0..2.0]) ∩ ProductDomain([0..1, 0..3]) isa Rectangle{Vector{Float64}}
        @test ProductDomain(0..1.0, 0..2.0) ∩ ProductDomain(0..1, 0..3) isa Rectangle{SVector{2,Float64}}

        # Generic functionality
        long_domain = ProductDomain([0..i for i in 1:20])
        show(io, long_domain)
        @test String(take!(io)) == "(0..1) × (0..2) × (0..3) × ... × (0..20)"
        @test isopenset(interior(UnitCube()))
        @test isclosedset(closure(interior(UnitCube())))
    end
    @testset "bounding box" begin
        @test boundingbox([0.2, -0.4, 1.0]) == -0.4..1.0
        @test boundingbox(Set([0.2, -0.4, 1.0])) == -0.4..1.0

        using DomainSets: unionbox, intersectbox

        @test unionbox(0..1) == 0..1
        @test unionbox(0..1, 2..3) == 0..3
        @test unionbox(0..1, 2..3, 5..6.0) === 0..6.0
        @test unionbox(ChebyshevInterval(), ChebyshevInterval()) === ChebyshevInterval()

        @test unionbox(EmptySpace{Float64}(), 2..3) === 2..3.0
        @test unionbox(2..3, EmptySpace{Float64}()) === 2..3.0
        @test unionbox(FullSpace{Int}(), 2..3.0) === FullSpace{Float64}()
        @test unionbox(2..3.0, FullSpace{Int}()) === FullSpace{Float64}()

        @test DomainSets.unionbox2(0..1, 0..1) == FullSpace()

        d1 = ProductDomain(0..1.0, 2..4.0)
        d2 = ProductDomain(0.5..1.5, 2.5..4.5)
        @test unionbox(d1, d2) == ProductDomain(0..1.5, 2..4.5)

        @test intersectbox(0..1) == 0..1
        @test intersectbox(0..1.5, 1..3) === 1.0..1.5
        @test intersectbox(0..1, 0.5..3, 0.0..4.0) === 0.5..1.0
        @test intersectbox(ChebyshevInterval(), ChebyshevInterval()) === ChebyshevInterval()

        @test intersectbox(EmptySpace{Int}(), 2..3.0) === EmptySpace{Float64}()
        @test intersectbox(2..3.0, EmptySpace{Int}()) === EmptySpace{Float64}()
        @test intersectbox(FullSpace{Float64}(), 2..3) === 2.0..3.0
        @test intersectbox(2..3, FullSpace{Float64}()) === 2.0..3.0

        @test DomainSets.intersectbox2(0..1, 0..1) == FullSpace()

        d1 = ProductDomain(0..1.0, 2..4.0)
        d2 = ProductDomain(0.5..1.5, 2.5..4.5)
        @test intersectbox(d1, d2) == ProductDomain(0.5..1.0, 2.5..4.0)

        @test boundingbox(MappedDomain(LinearMap(1/2), 2..3)) == 4.0..6.0
        @test boundingbox(MappedDomain(LinearMap(1/2), (2..3)^2)) == (4.0..6.0)^2
    end
end
