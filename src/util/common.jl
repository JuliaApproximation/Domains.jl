

isreal(::Type{<:Real}) = true
isreal(::Type{<:Complex}) = false
isreal(::Type{T}) where {T} = isreal(eltype(T))


#######################
# Composite structures
#######################

"""
Some types have composite structure, e.g. product domains, a union of domains.
These types contain a list of domains.

It is often undesirable to use `getindex` to access the elements of the composite
type. For this reason we introduce the `elements` functions. Composite types
can implement `elements` and provide a generic way to access their components.

`elements(t)`: returns the elements making up the composite type `t`

`element(t, i)`: return the `i`-th element of the composite type `t`

`numelements(t)`: return the number of elements of the composite type `t`
"""
elements(t) = ()

"""
Return the i-th element of a composite structure.

See also: `elements`.
"""
element(t, i) = elements(t)[i]

"""
Return the number of elements of a composite structure.

See also: `elements`.
"""
numelements(t) = length(elements(t))

"""
Is `x` composed of different elements?

See also: `elements`.
"""
iscomposite(t) = length(elements(t)) > 0

"Expand all arguments of type C into their components."
expand(::Type{C}) where {C} = ()
expand(::Type{C}, domain, domains...) where {C} =
    (domain, expand(C, domains...)...)
expand(::Type{C}, domain::C, domains...) where {C} =
    (elements(domain)..., expand(C, domains...)...)



#################
# Precision type
#################

"The floating point precision type associated with the argument."
prectype(x) = prectype(typeof(x))
prectype(::Type{<:Complex{T}}) where {T} = prectype(T)
prectype(::Type{<:AbstractArray{T}}) where {T} = prectype(T)
prectype(::Type{NTuple{N,T}}) where {N,T} = prectype(T)
prectype(::Type{Tuple{A}}) where {A} = prectype(A)
prectype(::Type{Tuple{A,B}}) where {A,B} = prectype(A,B)
prectype(::Type{Tuple{A,B,C}}) where {A,B,C} = prectype(A,B,C)
@generated function prectype(T::Type{<:Tuple{Vararg}})
    quote $(promote_type(map(prectype, T.parameters[1].parameters)...)) end
end
prectype(::Type{T}) where {T<:AbstractFloat} = T
prectype(::Type{T}) where {T<:Number} = prectype(float(T))

prectype(a...) = promote_type(map(prectype, a)...)

"Convert `x` such that its `prectype` equals `U`."
convert_prectype(x, ::Type{U}) where {U} = convert(to_prectype(typeof(x),U), x)

# function convert_prectype(x, ::Type{U}) where {U}
#     @warn "The order of arguments of convert_prectype has changed."
#     prectype(U, x)
# end

"Return the type to which `U` can be converted, such that the `prectype` becomes `T`."
to_prectype(::Type{T}, ::Type{U}) where {T,U} = error("Don't know how to convert the numtype of $(T) to $(U).")
to_prectype(::Type{T}, ::Type{U}) where {T <: Real,U <: Real} = U
to_prectype(::Type{Complex{T}}, ::Type{U}) where {T <: Real,U <: Real} = Complex{U}
to_prectype(::Type{SVector{N,T}}, ::Type{U}) where {N,T,U} = SVector{N,to_prectype(T,U)}
to_prectype(::Type{Vector{T}}, ::Type{U}) where {T,U} = Vector{to_prectype(T,U)}

"Promote the precision types of the arguments to a joined supertype."
promote_prectype(a) = a
promote_prectype(a, b) = _promote_prectype(prectype(a,b), a, b)
promote_prectype(a, b, c...) = _promote_prectype(prectype(a,b,c...), a, b, c...)
_promote_prectype(U, a) = convert_prectype(a, U)
_promote_prectype(U, a, b) = convert_prectype(a, U), convert_prectype(b, U)
_promote_prectype(U, a, b, c...) =
    (convert_prectype(a, U), convert_prectype(b, U), _promote_prectype(U, c...)...)

#################
# Numeric type
#################

"The numeric element type of x in a Euclidean space."
numtype(x) = numtype(typeof(x))
numtype(::Type{T}) where {T<:Number} = T
numtype(::Type{T}) where {T} = eltype(T)
numtype(::Type{NTuple{N,T}}) where {N,T} = T
numtype(::Type{Tuple{A,B}}) where {A,B} = promote_type(numtype(A), numtype(B))
numtype(::Type{Tuple{A,B,C}}) where {A,B,C} = promote_type(numtype(A), numtype(B), numtype(C))
numtype(::Type{Tuple{A,B,C,D}}) where {A,B,C,D} = promote_type(numtype(A), numtype(B), numtype(C), numtype(D))
@generated function numtype(T::Type{<:Tuple{Vararg}})
    quote $(promote_type(T.parameters[1].parameters...)) end
end

numtype(a...) = promote_type(map(numtype, a)...)

"Convert `x` such that its `numtype` equals `U`."
convert_numtype(x, ::Type{U}) where {U} = convert(to_numtype(typeof(x),U), x)

to_numtype(::Type{T}, ::Type{U}) where {T,U} = error("Don't know how to convert the numtype of $(T) to $(U).")
to_numtype(::Type{T}, ::Type{U}) where {T <: Number,U <: Number} = U
to_numtype(::Type{SVector{N,T}}, ::Type{U}) where {N,T,U} = SVector{N,to_numtype(T,U)}
to_numtype(::Type{Vector{T}}, ::Type{U}) where {T,U} = Vector{to_numtype(T,U)}

"Promote the numeric types of the arguments to a joined supertype."
promote_numtype(a) = a
promote_numtype(a, b) = _promote_numtype(numtype(a,b), a, b)
promote_numtype(a, b, c...) = _promote_numtype(numtype(a,b,c...), a, b, c...)
_promote_numtype(U, a) = convert_numtype(a, U)
_promote_numtype(U, a, b) = convert_numtype(a, U), convert_numtype(b, U)
_promote_numtype(U, a, b, c...) =
    (convert_numtype(a, U), convert_numtype(b, U), _promote_numtype(U, c...)...)



## Conversion from nested vectors to flat vectors and back

"""
Convert a vector from a cartesian format to a nested tuple according to the
given dimensions.

For example:
`convert_fromcartesian([1,2,3,4,5], Val{(2,2,1)}()) -> ([1,2],[3,4],5)`
"""
@generated function convert_fromcartesian(x::AbstractVector, ::Val{DIM}) where {DIM}
	dimsum = [0; cumsum([d for d in DIM])]
	E = Expr(:tuple, [ (dimsum[i+1]-dimsum[i] > 1 ? Expr(:call, :SVector, [:(x[$j]) for j = dimsum[i]+1:dimsum[i+1]]...) : :(x[$(dimsum[i+1])])) for i in 1:length(DIM)]...)
	return quote $(E) end
end

"The inverse function of `convert_fromcartesian`."
@generated function convert_tocartesian(x, ::Val{DIM}) where {DIM}
    dimsum = [0; cumsum([d for d in DIM])]
    E = vcat([[:(x[$i][$j]) for j in 1:DIM[i]] for i in 1:length(DIM)]...)
    quote SVector($(E...)) end
end
