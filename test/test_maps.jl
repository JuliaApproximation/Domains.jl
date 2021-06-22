
using DomainSets: ScalarAffineMap,
    VectorAffineMap,
    StaticAffineMap,
    GenericAffineMap,
    ScalarLinearMap,
    VectorLinearMap,
    StaticLinearMap,
    GenericLinearMap,
    ScalarTranslation,
    VectorTranslation,
    StaticTranslation,
    GenericTranslation,
    TupleProductMap


maps_to_test(T) = [
    StaticIdentityMap{T}(),
    VectorIdentityMap{T}(10),
    ConstantMap{T}(one(T)),
    ConstantMap{T}(SVector{2,T}(1,2)),
    ZeroMap{T}(),
    UnityMap{T}(),
    AffineMap(T(1.2), T(2.4)), # scalar map
    AffineMap(-T(1.2), T(2.4)), # scalar map with negative A
    AffineMap(randvec(T, 2, 2), randvec(T, 2)), # static map
    AffineMap(randvec(T, 3, 2), randvec(T, 3)), # static map, rectangular
    AffineMap(rand(T, 2, 2), rand(T, 2)), # vector map
    AffineMap(rand(T, 3, 2), rand(T, 3)), # vector map, rectangular
    AffineMap(LinearAlgebra.I, one(T)/2),  # use UniformScaling object as A
    AffineMap(LinearAlgebra.I, randvec(T,2)),  # use UniformScaling object as A
    DomainSets.GenericAffineMap(randvec(T, 2, 2), randvec(T, 2)),
    DomainSets.GenericAffineMap(T(1.2), randvec(T, 2)),
    DomainSets.GenericAffineMap(randvec(T, 3, 2), randvec(T, 3)),
    Translation(randvec(T, 3)),
    LinearMap(randvec(T, 2, 2)),
    LinearMap(randvec(T, 2)),
    LinearMap(randvec(T, 2, 2)) ∘ AffineMap(T(1.2), randvec(T, 2)),
    mapfrom_parameterdomain(T(5) * ComplexUnitCircle{T}() .+ T(2))
]

randvec(T,n) = SVector{n,T}(rand(n))
randvec(T,m,n) = SMatrix{m,n,T}(rand(m,n))

suitable_point_to_map(m::Map) = suitable_point_to_map(m, domaintype(m))

suitable_point_to_map(m::Map, ::Type{SVector{N,T}}) where {N,T} = SVector{N,T}(rand(N))
suitable_point_to_map(m::Map, ::Type{T}) where {T<:Number} = rand(T)
suitable_point_to_map(m::Map, ::Type{<:AbstractVector{T}}) where {T} = rand(T, mapsize(m,2))

suitable_point_to_map(m::DomainSets.ProductMap) =
    map(suitable_point_to_map, components(m))
suitable_point_to_map(m::DomainSets.VcatMap{N,T}) where {N,T} =
    SVector{N,T}(rand(T,N))

suitable_point_to_map(::CartToPolarMap{T}) where {T} = randvec(T,2)
suitable_point_to_map(::PolarToCartMap{T}) where {T} = randvec(T,2)

widertype(T) = widen(T)
widertype(::Type{SVector{N,T}}) where {N,T} = SVector{N,widen(T)}
widertype(::Type{Vector{T}}) where {T} = Vector{widen(T)}
widertype(::Type{Tuple{A}}) where {A} = Tuple{widen(A)}
widertype(::Type{Tuple{A,B}}) where {A,B} = Tuple{widen(A),widen(B)}
widertype(::Type{Tuple{A,B,C}}) where {A,B,C} = Tuple{widen(A),widen(B),widen(C)}
widertype(::Type{NTuple{N,T}}) where {N,T} = NTuple{N,widen(T)}

issquarematrix(A) = false
issquarematrix(A::AbstractArray) = size(A,1)==size(A,2)


function test_generic_inverse(m)
    M = mapsize(m,1)
    N = mapsize(m,2)
    x = suitable_point_to_map(m)
    y = m(x)

    if M==N
        inverse(m, y)      # trigger exception outside test if not implemented

        minv = inverse(m)
        @test minv(y) ≈ x
        @test inverse(m)(y) ≈ x
        @test inverse(m, y) ≈ x
        @test m\y ≈ x
        @test DomainSets.LazyInverse(m)(y) ≈ inverse(m, y)
    end
    if M >= N && numtype(m)!=BigFloat
        leftinverse(m, y)      # trigger exception outside test if not implemented

        mli = leftinverse(m)
        @test mli(y) ≈ x
        @test leftinverse(m, y) ≈ x
    end
    if M <= N && numtype(m)!=BigFloat
        rightinverse(m, y)      # trigger exception outside test if not implemented

        mri = rightinverse(m)
        @test m(mri(y)) ≈ y
        @test m(rightinverse(m, y)) ≈ y
    end
end

function test_generic_jacobian(m)
    x = suitable_point_to_map(m)
    # Trigger exception outside of test if jacobian(m, x) is not implemented.
    # test_generic_jacobian(m) should be called within a try/catch block
    jacobian(m, x)

    δ = sqrt(eps(prectype(m)))
    x2 = x .+ δ
    # we intentionally test jacobian(m, x) before testing jacobian(m)
    if !(m isa ProductMap)
        @test norm(m(x2) - m(x) + jacobian(m, x)*(x-x2)) < 100δ
    end
    jac = jacobian(m)
    @test jac(x) ≈ jacobian(m, x)
    j = jacobian(m, x)
    @test size(j) == mapsize(m)
    if j isa AbstractArray{Float64}
        y = similar(j)
        DomainSets.jacobian!(y, m, x)
        @test y ≈ jac(x)
    end
    if issquarematrix(jac(x))
        @test jacdet(m, x) ≈ det(jacobian(m, x))
    end
    if mapsize(m) == ()
        d = diffvolume(m, x)
        s = sqrt(det(transpose(jacobian(m,x))*jacobian(m,x)))
        @test d ≈ s || d ≈ -s
    else
        @test abs(diffvolume(m, x) - sqrt(det(transpose(jacobian(m,x))*jacobian(m,x)))) < 1e-7
    end
end

# generic map test suite
function test_generic_map(m)
    @test convert(Map{domaintype(m)}, m) == m

    x = suitable_point_to_map(m)
    @test applymap(m, x) == m(x)
    y = m(x)
    @test isreal(m) == (isreal(x) && isreal(y))

    S = domaintype(m)
    U = codomaintype(m)
    @test x isa S
    @test y isa U

    if isaffine(m) && !isconstant(m) && !(prectype(m) == BigFloat)
        test_generic_inverse(m)
    else
        try
            # The map may not support an inverse, let's try
            test_generic_inverse(m)
        catch e
            # make sure it failed because the inverse is not defined
            @test e isa MethodError || prectype(m) == BigFloat
        end
    end

    if isaffine(m)
        M = matrix(m)
        @test size(M) == mapsize(m)
        test_generic_jacobian(m)
    else
        try # jacobian may not be implemented
            test_generic_jacobian(m)
        catch e
        end
    end

    if domaintype(m) == Float64
        @test convert(Map{BigFloat}, m) isa Map{BigFloat}
        @test convert(Map{BigFloat}, m) == m
    end
    if prectype(m) == Float64
        U = widertype(domaintype(m))
        @test convert(Map{U}, m) isa Map{U}
        @test convert(Map{U}, m) == m
    end

    if isaffine(m)
        A = matrix(m)
        b = vector(m)
        @test m(x) == A*x+b
    end
end

function test_isomorphisms(T)
    m1 = DomainSets.VectorToNumber{T}()
    @test m1(SA[1.0]) == 1.0
    m1b = DomainSets.NumberToVector{T}()
    @test inverse(m1) == m1b
    @test inverse(m1, 0.4) == m1b(0.4)
    @test inverse(m1b) == m1
    @test inverse(m1b, SA[0.4]) == m1(SA[0.4])
    @test m1b(1.0) == SA[1.0]
    @test mapsize(m1) == (1,1)
    @test mapsize(m1b) == (1,)
    @test mapsize(m1 ∘ inverse(m1)) == ()
    @test mapsize(inverse(m1) ∘ m1) == (1,1)
    @test jacobian(m1 ∘ inverse(m1), 1) == one(T)
    @test jacobian(inverse(m1) ∘ m1, [1]) == ones(T,1,1)
    @test jacobian(m1) isa ConstantMap
    @test jacobian(m1b) isa ConstantMap

    m2 = DomainSets.VectorToComplex{T}()
    @test m2(SA[one(T), one(T)]) == 1 + im
    m2b = DomainSets.ComplexToVector{T}()
    @test inverse(m2) == m2b
    @test inverse(m2, zero(T)+0im) == m2b(zero(T)+0im)
    @test inverse(m2b) == m2
    @test inverse(m2b, SA[one(T),one(T)]) == one(T)+one(T)*im
    @test m2b(one(T)+one(T)*im) == SA[one(T),one(T)]
    @test mapsize(m2) == (1,2)
    @test mapsize(m2b) == (2,)
    @test mapsize(m2 ∘ m2b) == ()
    @test mapsize(m2b ∘ m2) == (2,2)

    m3 = DomainSets.VectorToTuple{2,T}()
    m3b = DomainSets.TupleToVector{2,T}()
    @test m3(SA[one(T), one(T)]) == (one(T),one(T))
    @test inverse(m3b, SA[one(T), one(T)]) == (one(T),one(T))
    @test inverse(m3) == m3b
    @test inverse(m3b) == m3
    @test m3b( (one(T),one(T)) ) == SA[one(T),one(T)]
    @test inverse(m3, (one(T),one(T)) ) == SA[one(T),one(T)]

    m4 = DomainSets.NestedToFlat{3,T,Tuple{Tuple{T,T},T},(2,1)}()
    m4b = DomainSets.FlatToNested{3,T,Tuple{Tuple{T,T},T},(2,1)}()
    x4 = ([T(1),T(2)],T(3))
    y4 = T[1,2,3]
    @test m4(x4) == y4
    @test m4b(y4) == x4
    @test inverse(m4) == m4b
    @test inverse(m4b) == m4
    @test inverse(m4, y4) == x4
    @test inverse(m4b, x4) == y4
end

function test_maps(T)
    generic_tests(T)

    # Test special maps
    test_identity_map(T)
    test_basic_maps(T)
    test_affine_maps(T)
    test_composite_map(T)
    test_product_map(T)
    test_wrapped_maps(T)
    test_scaling_maps(T)
    test_isomorphisms(T)
end

function generic_tests(T)
    for map in maps_to_test(T)
        @test prectype(map) == T
        test_generic_map(map)
    end
    # Test an additional composite map
    m1 = LinearMap(SMatrix{2,2}(1,2,3,4.0))
    m2 = CartToPolarMap()
    cmap = m1 ∘ m2 ∘ m1
    test_generic_map(cmap)
end

function test_affine_maps(T)
    A = rand(T,2,2)
    @test DomainSets.to_matrix(Vector{T}, A) == A
    @test DomainSets.to_matrix(T, 2) == 2
    @test DomainSets.to_matrix(SVector{2,T}, 2) == SMatrix{2,2}(2,0,0,2)
    @test DomainSets.to_matrix(SVector{2,T}, LinearAlgebra.I) == SMatrix{2,2}(1,0,0,1)
    @test DomainSets.to_matrix(Vector{T}, 2) == UniformScaling(2)
    @test DomainSets.to_matrix(Vector{T}, LinearAlgebra.I) == LinearAlgebra.I
    # test fallback with nonsensical call
    @test DomainSets.to_matrix(Tuple{Int}, 2) == 2

    @test DomainSets.to_matrix(T, A, 2) == A
    @test DomainSets.to_matrix(T, 2, 3) == 2
    @test DomainSets.to_matrix(T, UniformScaling(2), 3) == 2
    @test DomainSets.to_matrix(SVector{2,T}, 2, SVector(1,1)) == SMatrix{2,2}(2,0,0,2)
    @test DomainSets.to_matrix(Vector{T}, 2, [1,2]) == [2 0 ; 0 2]

    @test DomainSets.to_vector(T, 2) == 0
    @test DomainSets.to_vector(SVector{2,T}, 2) == SVector(0,0)
    @test DomainSets.to_vector(Vector{T}, A) == [0,0]
    @test DomainSets.to_vector(T, 2, 3) == 3

    if T != BigFloat    # BigFloat's make pinv fail for StaticArrays
        @test DomainSets.matrix_pinv(SMatrix{2,2}(rand(T),rand(T),rand(T),rand(T))) isa SMatrix{2,2}
        @test DomainSets.matrix_pinv(SVector(rand(T),rand(T))) isa Transpose{T,SVector{2,T}}
    end

    test_linearmap(T)
    test_translation(T)
    test_affinemap(T)
end



function test_linearmap(T)
    m1 = LinearMap(2one(T))
    @test m1 isa LinearMap{T}
    @test domaintype(m1) == T
    @test islinear(m1)
    @test isaffine(m1)
    @test m1(3) == 6
    @test m1(3one(T)) == 6
    @test m1(3one(widen(T))) isa widen(T)
    @test m1(SVector(1,2)) == SVector(2,4)
    @test m1([1.0,2.0]) == [2.0,4.0]
    mw1 = convert(Map{widen(T)}, m1)
    @test mw1 isa Map{widen(T)}
    @test mw1.A isa widen(T)
    @test jacobian(m1) isa ConstantMap{T}
    @test jacobian(m1, 1) == 2
    @test LinearMap(one(T)) == StaticIdentityMap{T}()

    m2 = LinearMap(2)
    @test domaintype(m2) == Int
    @test m2 === ScalarLinearMap(2)
    @test m2(one(T)) isa T
    @test m2 == convert(Map, 2)
    @test jacobian(m2, 1) == 2
    @test jacobian(m2) isa ConstantMap{Int}
    @test jacobian(m2, 1) == 2
    @test LinearMap(1) == StaticIdentityMap{T}()

    m3 = LinearMap(SMatrix{2,2}(one(T), 2one(T), 3one(T), 4one(T)))
    @test m3 isa LinearMap{SVector{2,T}}
    @test m3 === StaticLinearMap(m3.A)
    @test m3(SVector(1,2)) == SVector(7, 10)
    @test m3(SVector{2,T}(1,2)) == SVector{2,T}(7, 10)
    @test m3 ∘ m3 isa LinearMap
    @test LinearMap(SA[1 0; 0 1]) == StaticIdentityMap{domaintype(m3)}()

    A = rand(T,2,2)
    m4 = LinearMap(A)
    @test m4 isa LinearMap{Vector{T}}
    @test m4([1,2]) ==  A * [1,2]
    y = zeros(T,2)
    @test (DomainSets.applymap!(y, m4, [1,2]); y == A * [1,2])
    @test jacobian(m4, [1,2]) == A

    m5 = LinearMap{Vector{T}}(UniformScaling(2*one(T)))
    @test m5 isa LinearMap{Vector{T}}
    @test m5([1,2]) ==  2 * [1,2]
    @test jacobian(m5, [1,2]) == UniformScaling(2)

    m6 = LinearMap{SVector{2,T}}(2)
    @test m6 isa GenericLinearMap{SVector{2,T},T}
    @test m6(SA[one(T),one(T)]) == [2,2]
    @test leftinverse(m6) isa LinearMap{SVector{2,T}}
    @test leftinverse(m6)([2,2]) == [1,1]

    # Test construction and conversion
    @test LinearMap{T}(1) isa ScalarLinearMap{T}
    @test LinearMap{Vector{T}}(rand(Int,5,5)) isa VectorLinearMap{T}
    @test LinearMap{SVector{2,T}}(SMatrix{2,2}(1, 2, 3, 4)) isa StaticLinearMap{T,2,2,4}

    @test GenericLinearMap(one(T)) isa GenericLinearMap{T,T}
    @test GenericLinearMap(SMatrix{2,2,T}(1,1,1,1)) isa GenericLinearMap{SVector{2,T},SMatrix{2,2,T,4}}

    @test VectorLinearMap(rand(T,3,3)) isa VectorLinearMap{T}
    @test StaticLinearMap(SMatrix{2,2,T}(1,2,3,4)) isa StaticLinearMap{T,2,2,4}

    @test convert(Map{SVector{2,T}}, LinearMap(rand(T,2,2))) isa StaticLinearMap{T,2,2,4}
    @test convert(Map{T}, 1) isa LinearMap{T}
end



function test_translation(T)
    v = randvec(T,3)
    m = Translation(v)
    @test !islinear(m)
    @test isaffine(m)
    @test inverse(inverse(m)) == m
    @test jacobian(m) isa ConstantMap
    @test vector(jacobian(m)) == [1 0 0; 0 1 0; 0 0 1]
    @test jacdet(m) isa UnityMap{SVector{3,T},T}

    @test mapsize(Translation(one(T))) == ()

    # Test construction and conversion
    @test Translation(one(T)) isa ScalarTranslation{T}
    @test Translation{T}(1) isa ScalarTranslation{T}
    @test Translation(SVector{2,T}(1,2)) isa StaticTranslation{T,2}
    @test Translation{SVector{2,T}}(rand(T,2)) isa StaticTranslation{T,2}
    @test Translation(rand(Int,2)) isa VectorTranslation{Int}
    @test Translation{Vector{T}}(rand(Int,2)) isa VectorTranslation{T}
    @test Translation(1:10) isa GenericTranslation{Vector{Int}}
    @test Translation{Vector{T}}(1:10) isa GenericTranslation{Vector{T}}

    @test StaticTranslation(MVector(one(T),2)) isa StaticTranslation{T,2}
    @test StaticTranslation{T}(SVector(1,2)) isa StaticTranslation{T,2}

    @test VectorTranslation(1:10) isa VectorTranslation{Int}
    @test VectorTranslation{T}(1:10) isa VectorTranslation{T}

    @test GenericTranslation(1:10) isa GenericTranslation{Vector{Int}}
    @test GenericTranslation{Vector{T}}(1:10) isa GenericTranslation{Vector{T}}

    @test GenericTranslation(one(T)) isa GenericTranslation{T,T}
    @test GenericTranslation{T}(1) isa GenericTranslation{T,T}
    @test_throws InexactError GenericTranslation{Int}(one(T)/2)

    @test convert(Map{SVector{2,T}}, Translation(rand(T,2))) isa StaticTranslation{T,2}
end


function test_affinemap(T)
    m1 = AffineMap(T(2), T(3))
    @test m1 isa ScalarAffineMap{T}
    @test !islinear(m1)
    @test isaffine(m1)
    @test m1(2) == 7

    @test convert(ScalarAffineMap{BigFloat}, m1) == m1

    @test m1 ∘ m1 isa AffineMap
    @test (m1 ∘ m1)(2) == 2*(2*2+3)+3

    m2 = AffineMap(T(2), SVector{2,T}(1,2))
    @test m2 isa GenericAffineMap{SVector{2,T}}
    @test m2(SVector(1,2)) == SVector(3,6)
    @test mapsize(m2) == (2,2)

    m3 = AffineMap(UniformScaling(2*one(T)), [one(T),2*one(T)])
    @test m3 isa AffineMap{Vector{T}}
    @test mapsize(m3) == (2,2)
    @test m3([1,2]) ==  2 * [1,2] + [1,2]
    y = zeros(T,2)
    @test (DomainSets.applymap!(y, m3, [1,2]); y == m3([1,2]))
    @test jacobian(m3, [1,2]) == [2 0; 0 2]
    @test jacdet(m3, [1,2]) == 4

    @test AffineMap(LinearAlgebra.I, 2.0) isa ScalarAffineMap
    m4 = GenericAffineMap(LinearAlgebra.I, 2.0)
    @test m4.A isa UniformScaling{Float64}
    @test jacdet(m4, 3.0) == 1.0

    @test mapsize(AffineMap(rand(3),rand(3))) == (3,)
    @test mapsize(AffineMap(LinearAlgebra.I,2)) == ()

    # Test construction and conversion
    @test AffineMap(1, 2*one(T)) isa ScalarAffineMap{T}
    @test AffineMap{BigFloat}(1, 2.0) isa ScalarAffineMap{BigFloat}
    @test AffineMap(rand(T,2),rand(T,2)) isa GenericAffineMap{T}
    @test AffineMap(rand(T,2,2),rand(T,2)) isa VectorAffineMap{T}
    @test AffineMap{Vector{T}}(rand(Int,2,2),rand(Int,2)) isa VectorAffineMap{T}
    @test AffineMap{SVector{2,T}}(rand(T,2,2),rand(T,2)) isa StaticAffineMap{T,2,2,4}
    @test AffineMap{SVector{2,T}}(rand(Int,2,2),rand(T,2)) isa StaticAffineMap{T,2,2,4}
    @test AffineMap{SVector{2,T}}(rand(Int,2,2),rand(Int,2)) isa StaticAffineMap{T,2,2,4}
    @test AffineMap{SVector{2,T}}(1,rand(Int,2)) isa GenericAffineMap{SVector{2,T},T,SVector{2,T}}
    @test AffineMap{SVector{2,T}}(SMatrix{3,2}(3,2,1,4,5,6),SVector(1,2,3)) isa StaticAffineMap{T,2,3,6}

    @test GenericAffineMap(SMatrix{2,2,T}(1,2,3,4), [1,2]) isa GenericAffineMap{SVector{2,T},SMatrix{2,2,T,4}}
    @test GenericAffineMap(SMatrix{1,1,Int}(1),SVector{1,T}(1)) isa GenericAffineMap{SVector{1,T}}
    @test GenericAffineMap(rand(T,2,2),rand(Int,2)) isa GenericAffineMap{Vector{T},Matrix{T},Vector{T}}
    @test GenericAffineMap{Vector{T}}(1, [1,2]) isa GenericAffineMap{Vector{T},T,Vector{Int}}
    @test GenericAffineMap{Vector{T}}(rand(Int,2,2), rand(Int,2)) isa GenericAffineMap{Vector{T},Matrix{T},Vector{T}}
    @test convert(GenericAffineMap{Vector{T}}, GenericAffineMap(rand(Int,2,2), rand(Int,2))) isa GenericAffineMap{Vector{T}}
    @test GenericAffineMap(rand(T,2),rand(T,2)) isa GenericAffineMap{T}
    @test GenericAffineMap{T}(rand(Int,2),rand(Int,2)) isa GenericAffineMap{T}

    @test VectorAffineMap(rand(T,2,2),rand(T,2)) isa VectorAffineMap{T}

    @test StaticAffineMap(SMatrix{1,1,T}(1),SVector{1,T}(1)) isa StaticAffineMap{T}
    @test StaticAffineMap(SMatrix{1,1,T}(1),SVector{1,Int}(1)) isa StaticAffineMap{T}

    @test StaticAffineMap(rand(T,2,2),SVector{2,T}(1,1)) isa StaticAffineMap{T}
    @test_throws AssertionError StaticAffineMap(rand(T,3,2),SVector{2,T}(1,1))
    @test StaticAffineMap(SMatrix{2,2,T}(1,2,3,4),[1,1]) isa StaticAffineMap{T,2}
    @test StaticAffineMap(SMatrix{3,2,T}(1,2,3,4,5,6),SVector{3,T}(1,2,3)) isa StaticAffineMap{T,2,3,6}

    @test convert(Map{SVector{2,T}}, AffineMap(rand(2,2),rand(2))) isa StaticAffineMap{T,2,2,4}
end

function test_scaling_maps(T)
    test_generic_map(LinearMap(T(2)))
    test_generic_map(LinearMap(T(2), T(3)))
    test_generic_map(LinearMap(T(2), T(3), T(4)))
    test_generic_map(LinearMap(T(2), T(3), T(4), T(5)))
end

function test_identity_map(T)
    i1 = StaticIdentityMap{T}()
    i2 = StaticIdentityMap{SVector{2,T}}()
    test_generic_map(i1)
    test_generic_map(i2)
    @test i1 == i2
    @test hash(i1) == hash(i2)
    @test islinear(i1)
    @test isaffine(i1)
    @test convert(StaticIdentityMap{SVector{2,T}}, i1) === i2
    @test jacobian(i1) isa ConstantMap
    @test jacobian(i1, 1) == 1
    @test jacdet(i1, 1) == 1
    m1 = convert(DomainSets.AbstractAffineMap{T}, i1)
    @test m1 isa LinearMap{T}
    @test jacdet(m1, 1) == 1
    @test convert(DomainSets.AbstractAffineMap, i1) isa LinearMap{T}
    m2 = convert(DomainSets.LinearMap{T}, i1)
    @test m2 isa LinearMap{T}
    @test jacdet(m2, 1) == 1

    i3 = VectorIdentityMap{T}(10)
    test_generic_map(i3)
    r = rand(T, 10)
    @test i3(r) ≈ r
    @test hash(i3) == hash(VectorIdentityMap{Int}(10))

    @test IdentityMap() ∘ LinearMap(2) == LinearMap(2.0)
end

function test_basic_maps(T)
    @test UnityMap{SVector{2,Float64}}() == UnityMap{SVector{2,Float64},Float64}()
    @test hash(ConstantMap(2)) == hash(ConstantMap(2.0))
    @test DomainSets.absmap(ConstantMap(-2)) == ConstantMap(2)
end

function test_composite_map(T)
    a = T(0)
    b = T(1)
    c = T(2)
    d = T(3)
    ma = StaticIdentityMap{T}()
    mb = DomainSets.interval_map(a, b, c, d)

    r = suitable_point_to_map(ma)
    m1 = ma∘mb
    test_generic_map(m1)
    @test m1(r) ≈ ma(mb(r))
    m2 = m1∘mb
    test_generic_map(m2)
    @test m2(r) ≈ m1(mb(r))
    m3 = mb∘m2
    test_generic_map(m3)
    @test m3(r) ≈ mb(m2(r))
    m = m2∘m3
    test_generic_map(m)
    @test m(r) ≈ m2(m3(r))

    m5 = ComposedMap(LinearMap(rand(T,2,2)), AffineMap(rand(T,2,2),rand(T,2)))
    test_generic_map(m5)
    @test jacobian(m5) isa ConstantMap
    @test m5[Component(1)] isa LinearMap
    @test m5[Component(2)] isa AffineMap
    @test ComposedMap(m5[Component(1:2)]...) == m5
    @test_throws BoundsError m5[Component(3)]

    m6 = DomainSets.multiply_map(ma,ma)
    @test m6(one(T)/2) == one(T)/4
    @test jacobian(m6) isa DomainSets.SumMap
    @test jacobian(m6)(one(T)) == 2
    @test jacobian(m6, one(T)) == 2

    m7 = DomainSets.sum_map(ma,ma)
    @test m7(one(T)) == 2
    @test jacobian(m7) isa ConstantMap
    @test jacobian(m7, one(T)) == 2

    @test mapsize(ComposedMap(LinearMap(2),LinearMap(rand(T,2)),LinearMap(rand(T,2,2)))) == (2,)
    @test mapsize(ComposedMap(LinearMap(rand(T,2)),LinearMap(rand(T,2)'),LinearMap(rand(T,2)))) == (2,)
    @test mapsize(ComposedMap(LinearMap(rand(T,2,2)),LinearMap(rand(T,2)'),LinearMap(2))) == (1,2)
    @test mapsize(ComposedMap(LinearMap(one(T)),LinearMap(one(T)))) == ()

    @test DomainSets.composedmap() == ()
    @test DomainSets.composedmap(ma) == ma
    @test DomainSets.composedmap(ma,ma) == ma
    @test DomainSets.composedmap(ma,ma,ma) == ma

    @test DomainSets.composite_jacobian(ma) == jacobian(ma)

    @test DomainSets.multiply_map() == ()
    @test DomainSets.multiply_map(ma) == ma
    @test DomainSets.multiply_map(ma,ma)(2*one(T)) == 4
    @test DomainSets.multiply_map(ma,ma,ma)(2*one(T)) == 8

    @test DomainSets.sum_jacobian() == ()
    @test DomainSets.sum_jacobian(ma) == jacobian(ma)
end

function test_product_map(T)
    ma = StaticIdentityMap{T}()
    mb = DomainSets.interval_map(T(0), T(1), T(2), T(3))

    r1 = suitable_point_to_map(ma)
    r2 = suitable_point_to_map(ma)
    r3 = suitable_point_to_map(ma)
    r4 = suitable_point_to_map(ma)
    r5 = suitable_point_to_map(ma)

    m1 = productmap(ma,mb)
    test_generic_map(m1)
    @test m1(SVector(r1,r2)) ≈ SVector(ma(r1),mb(r2))
    m2 = productmap(m1,mb)
    test_generic_map(m2)
    @test m2(SVector(r1,r2,r3)) ≈ SVector(ma(r1),mb(r2),mb(r3))
    m3 = productmap(mb,m2)
    test_generic_map(m3)
    @test m3(SVector(r1,r2,r3,r4)) ≈ SVector(mb(r1),ma(r2),mb(r3),mb(r4))
    m4 = productmap(m1,m2)
    test_generic_map(m4)
    @test m4(SVector(r1,r2,r3,r4,r5)) ≈ SVector(m1(SVector(r1,r2))...,m2(SVector(r3,r4,r5))...)

    @test ProductMap(ma, mb) == ProductMap([ma,mb])
    @test hash(ProductMap(ma, mb)) == hash(ProductMap([ma,mb]))

    m5 = productmap(AffineMap(SMatrix{2,2,T}(1.0,2,3,4), SVector{2,T}(1,3)), LinearMap{T}(2.0))
    @test domaintype(component(m5,1)) == SVector{2,T}
    @test domaintype(component(m5,2)) == T
    test_generic_map(m5)
    x = SVector{3,T}(rand(T,3))
    @test m5(x) ≈ SVector(component(m5,1)(SVector(x[1],x[2]))...,component(m5,2)(x[3]))

    m6 = ProductMap([ma,mb])
    @test m6 isa DomainSets.VectorProductMap
    @test mapsize(m6) == (2,2)
    @test convert(Map{SVector{2,T}}, m6) isa DomainSets.VcatMap
    test_generic_map(m6)

    d1 = ProductDomain(T(1)..T(2), T(1)..T(2))
    d2 = ProductDomain(T(2)..T(4), T(2)..T(4))
    @test mapto(d1, d2) isa DomainSets.VcatMap
    @test jacobian(mapto(d1,d2)) isa ConstantMap
    @test jacdet(mapto(d1,d2)) == ConstantMap{SVector{2,T}}(4)
    d1v = ProductDomain([T(1)..T(2), T(1)..T(2)])
    d2v = ProductDomain([T(2)..T(4), T(2)..T(4)])
    @test mapto(d1v, d2v) isa DomainSets.VectorProductMap
    @test jacobian(mapto(d1v,d2v)) isa ConstantMap
    @test jacdet(mapto(d1v,d2v)) == ConstantMap{Vector{T}}(4)

    @test components(productmap(LinearMap(1), LinearMap(one(T)))) isa Tuple{<:Map{T},<:Map{T}}
end

using DomainSets: WrappedMap
function test_wrapped_maps(T)
    m1 = WrappedMap{T}(cos)
    m2 = WrappedMap{T}(sin)
    @test m1(one(T)) ≈ cos(one(T))
    @test m2(one(T)) ≈ sin(one(T))
    @test m1 == cos
    @test sin == m2
    m3 = m1 ∘ m2
    @test m3(one(T)) ≈ cos(sin(one(T)))

    @test WrappedMap(cos) isa WrappedMap{Float64}
    @test convert(Map, cos) isa WrappedMap
    @test convert(Map{BigFloat}, m1) isa WrappedMap{BigFloat}

    @test convert(Map{T}, cos) isa DomainSets.WrappedMap{T,typeof(cos)}

    @test convert(Map, LinearAlgebra.I) isa GenericLinearMap{Vector{Any}}
    @test convert(Map{Vector{T}}, LinearAlgebra.I) isa GenericLinearMap{Vector{T}}
end

function test_mixed_maps()
    m1 = composedmap(cos, sin)
    @test domaintype(m1) == Any
    @test m1(0.4) == sin(cos(0.4))

    m2 = composedmap(AffineMap(2.0, 3.0), cos)
    @test domaintype(m2) == Float64
    @test m2(0.4) ≈ cos(2*0.4+3)
    @inferred m2(0.4)
    @test repr(m2) == "cos ∘ (x -> 2.0 * x + 3.0)"

    m3 = composedmap(sin, AffineMap(2.0, 3.0), cos)
    @test m3(0.5) ≈ cos(2*sin(0.5)+3)
    @test repr(m3) == "cos ∘ (x -> 2.0 * x + 3.0) ∘ sin"
    @test domaintype(m3) == Any
    @inferred m3(0.5)

    m4 = productmap(sin, cos)
    @test m4 isa TupleProductMap
    @test domaintype(m4) == Tuple{Any,Any}
    @test m4(0.3,0.5) == (sin(0.3), cos(0.5))
end

@testset "maps" begin
    test_maps(Float64)
    test_maps(BigFloat)
    test_mixed_maps()
end
