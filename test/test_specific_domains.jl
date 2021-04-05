using StaticArrays, DomainSets, Test

import DomainSets: MappedDomain, similar_interval

struct Basis3Vector <: StaticVector{3,Float64} end

Base.getindex(::Basis3Vector, k::Int) = k == 1 ? 1.0 : 0.0

const io = IOBuffer()


struct NamedHyperBall <: DomainSets.DerivedDomain{SVector{2,Float64}}
    domain  ::  Domain{SVector{2,Float64}}

    NamedHyperBall() = new(2UnitDisk())
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
        d2 = 0..1
        @test d1 ∩ d2 == d1
        @test d2 ∩ d1 == d1
        @test d1 ∪ d2 == d2
        @test d2 ∪ d1 == d2
        @test d1 \ d2 == d1
        @test d2 \ d1 == d2
        @test d2 \ d2 == d1

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
        d2 = 0..1
        @test d1 ∪ d2 == d1
        @test d1 ∩ d2 == d2
        @test d2 ∩ d1 == d2
        @test typeof(FullSpace(0..1) .+ 1) <: FullSpace
        @test typeof(FullSpace(0..1) * 3) <: FullSpace
        @test infimum(d1) == typemin(Float64)
        @test supremum(d1) == typemax(Float64)

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
        @test isclosedset(d)
        @test !isopenset(d)
        @test dimension(d) == 1
        @test isempty(interior(d))
        @test closure(d) == d

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

        @test dimension(Point([1,2,3]))==3
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

            @test similar_interval(0..1, 0, big(1.0)) isa ClosedInterval{BigFloat}
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
        end
        @testset "ChebyshevInterval{$T}" begin
            d = ChebyshevInterval{T}()
            @test leftendpoint(d) == -one(T)
            @test rightendpoint(d) == one(T)
            @test minimum(d) == infimum(d) == leftendpoint(d)
            @test maximum(d) == supremum(d) == rightendpoint(d)

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

    @testset "unit ball" begin
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

        @test convert(SublevelSet, UnitDisk()) isa SublevelSet{SVector{2,Float64},:closed}
        @test convert(SublevelSet, EuclideanUnitBall{2,Float64,:open}()) isa SublevelSet{SVector{2,Float64},:open}

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
        @test SA[1.4, 1.4] ∈ D
        @test SA[1.5, 1.5] ∉ D
        @test typeof(1.2 * D)==typeof(D * 1.2)
        @test SA[1.5,1.5] ∈ 1.2 * D
        @test SA[1.5,1.5] ∈ D * 1.2
        @test !isempty(D)
        # TODO: implement and test isclosedset and isopenset for mapped domains

        D = 2UnitDisk() .+ SA[1.0,1.0]
        @test SA[2.4, 2.4] ∈ D
        @test SA[3.5, 2.5] ∉ D
        @test !isempty(D)

        B = UnitBall()
        @test SA[1.,0.0,0.] ∈ B
        @test SA[1.,0.1,0.] ∉ B
        @test !isempty(B)
        @test isclosedset(B)
        @test !isopenset(B)
        @test boundary(B) == UnitSphere()
        @test isopenset(interior(B))

        B = 2UnitBall()
        @test SA[1.9,0.0,0.0] ∈ B
        @test SA[0,-1.9,0.0] ∈ B
        @test SA[0.0,0.0,-1.9] ∈ B
        @test SA[1.9,1.9,0.0] ∉ B
        @test !isempty(B)

        B = 2.0UnitBall() .+ SA[1.0,1.0,1.0]
        @test SA[2.9,1.0,1.0] ∈ B
        @test SA[1.0,-0.9,1.0] ∈ B
        @test SA[1.0,1.0,-0.9] ∈ B
        @test SA[2.9,2.9,1.0] ∉ B
        @test !isempty(B)

        C = VectorUnitBall(4)
        @test [1, 0, 0, 0] ∈ C
        @test [0.0,0.1,0.2,0.1] ∈ C
        @test SA[0.0,0.1,0.2,0.1] ∈ C
        @test_logs (:warn, "`in`: incompatible combination of vector with length 2 and domain 'the 4-dimensional closed unit ball' with dimension 4. Returning false.") [0.0,0.1] ∉ C
        @test [0.0,1.1,0.2,0.1] ∉ C
        @test !isempty(C)
        @test isclosedset(C)
        @test !isopenset(C)
        @test boundary(C) == VectorUnitSphere(4)
        @test dimension(C) == 4
        @test isopenset(interior(C))

        D = VectorUnitBall{Float64,:open}(4)
        @test !in([1, 0, 0, 0], D)
        @test in([1-eps(Float64), 0, 0, 0], D)
        @test approx_in([1.1, 0, 0, 0], D, 0.2)
        @test !approx_in([1.1, 0, 0, 0], D, 0.01)
        @test !isempty(D)
        @test approx_in(SA[1.01,0.0,0.0,0.0], D, 0.05)

        @test isclosedset(DomainSets.StaticUnitBall{SVector{2,Float64}}())
        @test EuclideanUnitBall{2}() isa EuclideanUnitBall{2,Float64}

        show(io, EuclideanUnitBall{2,Float64,:open}())
        @test String(take!(io)) == "the 2-dimensional open unit ball"
        show(io, UnitCircle())
        @test String(take!(io)) == "the unit circle"
    end

    @testset "custom named ball" begin
        B = NamedHyperBall()
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

        show(io,C)
        @test String(take!(io)) == "the complex unit circle (T=Complex{Float64})"
        show(io,D)
        @test String(take!(io)) == "the complex unit disk (T=Complex{Float64})"
        show(io,D2)
        @test String(take!(io)) == "the complex open unit disk (T=Complex{BigFloat})"

        @test pseudolevel(ComplexUnitCircle(), 0.1) isa SublevelSet{Complex{Float64},:open}
        p = pseudolevel(ComplexUnitCircle(), 0.1)
        @test 0.8 ∉ p
        @test 0.95 ∈ p
        @test 1+0.1im ∈ p
        @test 1.1+0.2im ∉ p
    end

    @testset "sphere" begin
        C = UnitCircle()
        @test SA[1.,0.] ∈ C
        @test SA[1.,1.] ∉ C
        @test approx_in(SA[1.,0.], C)
        @test !approx_in(SA[1.,1.], C)
        @test !isempty(C)
        @test isclosedset(C)
        @test !isopenset(C)
        p = parameterization(C)
        x = applymap(p, 1/2)
        @test gradient(p, 0.4) ≈ SA[-2pi*sin(2pi*0.4), 2pi*cos(2pi*0.4)]
        @test approx_in(x, C)
        q = leftinverse(p)
        @test applymap(q, x) ≈ 1/2
        @test applymap(q, -x) ≈ 1
        @test rightinverse(q) == p

        @test convert(LevelSet, UnitCircle()) isa LevelSet{SVector{2,Float64}}
        @test convert(LevelSet{SVector{2,BigFloat}}, UnitCircle()) isa LevelSet{SVector{2,BigFloat}}
        @test pseudolevel(UnitCircle(), 0.1) isa SublevelSet
        @test SA[1.05,0] ∈ pseudolevel(UnitCircle(), 0.1)
        @test SA[1.15,0] ∉ pseudolevel(UnitCircle(), 0.1)

        C2 = convert(Domain{SVector{2,BigFloat}}, C)
        @test eltype(C2) == SVector{2,BigFloat}

        C = 2UnitCircle() .+ SA[1.,1.]
        @test approx_in(SA[3.,1.], C)

        C = UnitCircle() .+ SA[1,1]
        @test approx_in(SA[2,1], C)

        S = UnitSphere()
        @test SA[1.,0.,0.] ∈ S
        @test SA[1.,0.,1.] ∉ S
        @test approx_in(SA[cos(1.),sin(1.),0.], S)
        @test !isempty(S)
        S2 = convert(Domain{SVector{3,BigFloat}}, S)
        @test eltype(S2) == SVector{3,BigFloat}

        @test Basis3Vector() in S

        S = 2 * UnitSphere() .+ SA[1.,1.,1.]
        @test approx_in(SA[1. + 2*cos(1.),1. + 2*sin(1.),1.], S)
        @test !approx_in(SA[4.,1.,5.], S)

        D = UnitCircle()
        @test convert(Domain{SVector{2,BigFloat}}, D) ≡ UnitCircle{BigFloat}()
        @test SVector(1,0) in D
        @test SVector(nextfloat(1.0),0) ∉ D

        D = UnitSphere()
        @test convert(Domain{SVector{3,BigFloat}}, D) ≡ UnitSphere{BigFloat}()
        @test SVector(1,0,0) in D
        @test SVector(nextfloat(1.0),0,0) ∉ D
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

        C = DomainSets.cylinder()
        @test eltype(C) == SVector{3,Float64}
        C2 = DomainSets.cylinder(1.0, 2)
        @test SA[0.5,0.2,1.5] ∈ C2
    end

    @testset "mapped_domain" begin
        @test MappedDomain(0..1.0, cos) isa MappedDomain{Float64}
        @test MappedDomain{Float64}(0..1.0, cos) isa MappedDomain{Float64}
        # Test chaining of maps
        D = UnitCircle()
        D1 = 2 * D
        @test typeof(D1) <: MappedDomain
        @test typeof(superdomain(D1)) <: UnitHyperSphere
        @test isclosedset(D1)
        @test !isopenset(D1)
        @test convert(Domain{SVector{2,BigFloat}}, D1) isa MappedDomain{SVector{2,BigFloat}}
        D2 = 2 * D1
        @test typeof(superdomain(D2)) <: UnitHyperSphere

        D = UnitInterval()^2
        show(io,rotate(D,1.))
        @test String(take!(io)) == "A mapped domain based on 0.0..1.0 (Unit) x 0.0..1.0 (Unit)"

        D = rotate(UnitInterval()^2, π)
        @test SA[-0.9, -0.9] ∈ D
        @test SA[-1.1, -1.1] ∉ D

        D = rotate(UnitInterval()^2, π, SA[-.5,-.5])
        @test SA[-1.5, -1.5] ∈ D
        @test SA[-0.5, -0.5] ∉ D

        D = rotate(UnitInterval()^3 .+ SA[-.5,-.5,-.5], pi, pi, pi)
        @test SA[0.4, 0.4, 0.4] ∈ D
        @test SA[0.6, 0.6, 0.6] ∉ D

        D = rotate((-1.5.. 2.2) × (0.5 .. 0.7) × (-3.0 .. -1.0), π, π, π, SA[.35, .65, -2.])
        @test SA[0.9, 0.6, -2.5] ∈ D
        @test SA[0.0, 0.6, 0.0] ∉ D

        B = 2 * VectorUnitBall(10)
        @test dimension(B) == 10
    end

    @testset "simplex" begin
        d = UnitSimplex{2}()
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

        @test approx_in(SA[-0.1,-0.1], d, 0.1)
        @test !approx_in(SA[-0.1,-0.1], d, 0.09)

        @test corners(d) == [ SA[0.0,0.0], SA[1.0,0.0], SA[0.0,1.0]]

        @test convert(Domain{SVector{2,BigFloat}}, d) == UnitSimplex{2,BigFloat}()

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
        @test SA[0.1,0.3] ∉ d4
        @test SA[-0.1,0.3] ∈ d4
        @test SA[-0.1,-0.3] ∉ d4
        @test SA[-0.1,0.3] ∈ d4_open
        convert(Domain{SVector{2,BigFloat}}, d4) isa SubzeroSet{SVector{2,BigFloat},:closed}
        convert(Domain{SVector{2,BigFloat}}, d4_open) isa SubzeroSet{SVector{2,BigFloat},:open}

        d5 = SuperlevelSet(cos, 0.5)
        d5_open = SuperlevelSet{Float64,:open}(cos, 0.5)
        @test d5 isa SuperlevelSet{Float64,:closed}
        @test interior(d5) == d5_open
        @test closure(d5_open) == d5
        @test closure(d5) == d5
        @test interior(d5_open) == d5_open
        @test 3.0 ∉ d5
        @test 0.0 ∈ d5
        @test 0.0 ∈ d5
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

        @test convert(IndicatorFunction, 0..1) isa IndicatorFunction
        @test convert(IndicatorFunction, d) == d
        @test convert(Domain{BigFloat}, d) isa IndicatorFunction{BigFloat}
        @test 0.5 ∈ convert(IndicatorFunction, 0..1)

        d2 = DomainSets.UntypedIndicatorFunction(ispositive)
        @test d2 isa DomainSets.UntypedIndicatorFunction{Float64}
        @test DomainSets.indicatorfunction(d2) == ispositive
        @test convert(Domain{BigFloat}, d2) isa DomainSets.UntypedIndicatorFunction{BigFloat}
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
        @test convert(Domain{SVector{2,BigFloat}}, d1) isa VcatDomain
        d1w = convert(Domain{SVector{2,BigFloat}}, d1)
        @test eltype(d1w) == SVector{2,BigFloat}
        @test eltype(element(d1w, 1)) == BigFloat
        @test eltype(element(d1w, 2)) == BigFloat

        @test VcatDomain( (-1..1, -2..2)) isa VcatDomain{2,Int,(1,1),Tuple{ClosedInterval{Int64}, ClosedInterval{Int64}}}

        show(io,d1)
        @test String(take!(io)) == "-1.0..1.0 x -1.0..1.0"

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
        @test_logs (:warn, "`in`: incompatible combination of vector with length 3 and domain '-1.0..1.0 x -1.0..1.0' with dimension 2. Returning false.") SA[0.0,0.0,0.0] ∉ d1
        @test_logs (:warn, "`in`: incompatible combination of vector with length 1 and domain '-1.0..1.0 x -1.0..1.0' with dimension 2. Returning false.") SA[0.0] ∉ d1
        @test_logs (:warn, "`in`: incompatible combination of vector with length 3 and domain '-1.0..1.0 x -1.0..1.0' with dimension 2. Returning false.") [0.0,0.0,0.0] ∉ d1
        @test_logs (:warn, "`in`: incompatible combination of vector with length 1 and domain '-1.0..1.0 x -1.0..1.0' with dimension 2. Returning false.") [0.0] ∉ d1

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
        @test eltype(element(d,1)) == Float64
        @test eltype(element(d,2)) == Float64
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
        @test eltype(element(d1big,1)) == BigFloat
        @test eltype(element(d1big,2)) == BigFloat

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
        @test UnitHyperCube(Val(2)) isa EuclideanUnitCube{2,Float64}
        @test UnitHyperCube{BigFloat}(Val(2)) isa EuclideanUnitCube{2,BigFloat}

        d1 = VectorUnitCube{Float64}(4)
        @test VectorUnitCube(4) == d1
        @test UnitHyperCube(4) == d1
        @test UnitHyperCube{Float64}(4) == d1
        @test dimension(d1) == 4
        @test element(d1, 1) == 0..1
        @test SA[0.9,0.9,0.4,0.2] ∈ d1
        @test [1.2,0.3,0.4,0.6] ∉ d1

        @test ProductDomain([UnitInterval(),UnitInterval()]) isa DomainSets.DynamicUnitCube{Float64}
        @test ProductDomain{Vector{BigFloat}}([UnitInterval(),UnitInterval()]) isa VectorProductDomain{Vector{BigFloat},Vector{UnitInterval{BigFloat}}}
        @test ProductDomain{SVector{2,BigFloat}}(UnitInterval(),UnitInterval()) isa EuclideanUnitCube{2,BigFloat}
        @test ProductDomain{SVector{2,BigFloat}}(SVector(UnitInterval(),UnitInterval())) isa EuclideanUnitCube{2,BigFloat}

        D = UnitInterval()^2
        @test D isa EuclideanUnitCube{2,Float64}
        @test SA[0.9, 0.9] ∈ D
        @test SA[1.1, 1.1] ∉ D
        @test !isempty(D)
        @test isclosedset(D)
        @test !isopenset(D)

        @test approx_in(SA[-0.1,-0.1], D, 0.1)
        @test !approx_in(SA[-0.1,-0.1], D, 0.09)

        #Cube
        D = (-1.5 .. 2.2) × (0.5 .. 0.7) × (-3.0 .. -1.0)
        @test SA[0.9, 0.6, -2.5] ∈ D
        @test SA[0.0, 0.6, 0.0] ∉ D
    end
    @testset "HyperRectangle" begin
        d1 = (-1.0..1.0) × (-1.0..1.0)

        d4 = d1 × (-1.0..1.0)
        @test d4 isa HyperRectangle
        @test SA[0.5,0.5,0.8] ∈ d4
        @test SA[-1.1,0.3,0.1] ∉ d4
        @test point_in_domain(d4) ∈ d4

        d5 = (-1.0..1.)×d1
        @test d5 isa HyperRectangle
        @test SA[0.,0.5,0.5] ∈ d5
        @test SA[0.,-1.1,0.3] ∉ d5
        @test point_in_domain(d5) ∈ d5

        d6 = d1 × d1
        @test d6 isa HyperRectangle
        @test SA[0.,0.,0.5,0.5] ∈ d6
        @test SA[0.,0.,-1.1,0.3] ∉ d6
        @test point_in_domain(d6) ∈ d6

        @test HyperRectangle( SA[1,2], SA[2.0,3.0]) isa HyperRectangle{SVector{2,Float64}}
        @test HyperRectangle([0..1, 2..3]) isa HyperRectangle{Vector{Int}}

        @test_throws ErrorException HyperRectangle(UnitCircle(), UnitDisk())
        @test_throws ErrorException HyperRectangle(OpenInterval(1,2), 3..4)

        bnd = boundary(HyperRectangle([1,2],[3,4]))
        @test [1,3] ∈ bnd
        @test [1,2.5] ∈ bnd
        @test [1.5,4] ∈ bnd
        @test [1.5,3.5] ∉ bnd
    end
    @testset "fixed product domains" begin
        d1 = ProductDomain(ChebyshevInterval(), ChebyshevInterval())
        @test d1 isa DomainSets.FixedIntervalProduct
        @test d1 isa DomainSets.ChebyshevProductDomain
        @test element(d1,1) isa ChebyshevInterval{Float64}
        @test elements(d1) isa NTuple{2,ChebyshevInterval{Float64}}
        @test convert(Domain{SVector{2,BigFloat}}, d1) isa DomainSets.ChebyshevProductDomain{2,BigFloat}
        @test ProductDomain(elements(d1)) == d1
        @test ProductDomain{SVector{2,BigFloat}}(elements(d1)) isa DomainSets.ChebyshevProductDomain{2,BigFloat}
        @test ProductDomain{SVector{2,BigFloat}}(elements(d1)...) isa DomainSets.ChebyshevProductDomain{2,BigFloat}
    end
    @testset "ProductDomain" begin
        d1 = 0..1.0
        d2 = 0..2.0
        d3 = UnitCircle()
        @test ProductDomain{SVector{2,Float64}}(d1, d2) isa VcatDomain
        @test ProductDomain{Tuple{Float64,Float64}}(d1, d2) isa TupleProductDomain
        @test ProductDomain{Vector{Float64}}([d1; d2]) isa VectorProductDomain

        @test ProductDomain((d1,d3)) isa VcatDomain
        @test ProductDomain((d1,d2)) isa HyperRectangle

        @test volume(ProductDomain(d1,d2)) == 2

        @test ProductDomain(SVector(0..1, 0..2)) isa HyperRectangle{SVector{2,Int}}
        @test ProductDomain(1.05 * UnitDisk(), -1.0 .. 1.0) isa VcatDomain{3,Float64}
        @test ProductDomain(['a','b'], 0..1) isa TupleProductDomain

        @test ProductDomain{Tuple{Float64,Float64}}(0..0.5, 0..0.7) isa TupleProductDomain{Tuple{Float64,Float64}}
        # Some conversions
        @test convert(Domain{Vector{Float64}}, (-1..1)^2) isa VectorProductDomain{Vector{Float64}}
        @test convert(Domain{SVector{2,Float64}}, ProductDomain([-1..1,-2..2])) isa VcatDomain{2,Float64}
        @test convert(Domain{Vector{Float64}}, TupleProductDomain(-1..1,-2..2)) isa VectorDomain{Float64}

        # intersection of product domains
        @test ProductDomain([0..1.0, 0..2.0]) ∩ ProductDomain([0..1, 0..3]) isa HyperRectangle{Vector{Float64}}
        @test ProductDomain(0..1.0, 0..2.0) ∩ ProductDomain(0..1, 0..3) isa HyperRectangle{SVector{2,Float64}}

        # Generic functionality
        long_domain = ProductDomain([0..i for i in 1:20])
        show(io, long_domain)
        @test String(take!(io)) == "0..1 x 0..2 x 0..3 x 0..4 x 0..5 x ... x 0..16 x 0..17 x 0..18 x 0..19 x 0..20"
        @test isopenset(interior(UnitCube()))
        @test isclosedset(closure(interior(UnitCube())))
    end
end
