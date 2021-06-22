
"A `FunctionLevelSet` is a set that derives from the levels of a function."
abstract type FunctionLevelSet{T} <: Domain{T} end

# Convenience: assume the function and the level are stored fields
levelfun(d::FunctionLevelSet) = d.f
levelfun(d::FunctionLevelSet, x) = d.f(x)
level(d::FunctionLevelSet) = d.level


"Supertype of level set domains of the form `f(x)=C`."
abstract type AbstractLevelSet{T} <: FunctionLevelSet{T} end

indomain(x, d::AbstractLevelSet) = levelfun(d, x) == level(d)

show(io::IO, d::AbstractLevelSet) =
    print(io, "level set f(x) = $(level(d)) with f = $(levelfun(d))")

==(d1::AbstractLevelSet, d2::AbstractLevelSet) = levelfun(d1)==levelfun(d2) &&
    level(d1)==level(d2)
hash(d::AbstractLevelSet, h::UInt) = hashrec("AbstractLevelSet", levelfun(d), level(d), h)

"The domain defined by `f(x)=0` for a given function `f`."
struct ZeroSet{T,F} <: AbstractLevelSet{T}
    f   ::  F
end

ZeroSet(f) = ZeroSet{Float64}(f)
ZeroSet{T}(f::F) where {T,F} = ZeroSet{T,F}(f)

level(d::ZeroSet) = 0

similardomain(d::ZeroSet, ::Type{T}) where {T} = ZeroSet{T}(levelfun(d))


"The domain defined by `f(x)=C` for a given function `f` and constant `C`."
struct LevelSet{T,F,S} <: AbstractLevelSet{T}
    f       ::  F
    level   ::  S
end

LevelSet(f, level) = LevelSet{typeof(level)}(f, level)
LevelSet{T}(f::F, level::S) where {T,F,S} = LevelSet{T,F,S}(f, level)

similardomain(d::LevelSet, ::Type{T}) where {T} = LevelSet{T}(levelfun(d), level(d))

convert(::Type{LevelSet}, d::ZeroSet{T}) where {T} = LevelSet{T}(levelfun(d), level(d))
convert(::Type{LevelSet{T}}, d::ZeroSet) where {T} = LevelSet{T}(levelfun(d), level(d))


"Supertype of sublevel set domains."
abstract type AbstractSublevelSet{T,C} <: FunctionLevelSet{T} end

indomain(x, d::AbstractSublevelSet{T,:closed}) where {T} = levelfun(d, x) <= level(d)
indomain(x, d::AbstractSublevelSet{T,:open}) where {T} = levelfun(d, x) < level(d)

show(io::IO, d::AbstractSublevelSet{T,:closed}) where {T} =
    print(io, "sublevel set f(x) <= $(level(d)) with f = $(levelfun(d))")
show(io::IO, d::AbstractSublevelSet{T,:open}) where {T} =
    print(io, "sublevel set f(x) < $(level(d)) with f = $(levelfun(d))")

==(d1::AbstractSublevelSet, d2::AbstractSublevelSet) = levelfun(d1)==levelfun(d2) &&
    level(d1)==level(d2)
hash(d::AbstractSublevelSet, h::UInt) =
    hashrec("AbstractSublevelSet", levelfun(d), level(d), h)

"The domain where `f(x) <= 0` (or `f(x) < 0`)."
struct SubzeroSet{T,C,F} <: AbstractSublevelSet{T,C}
    f   ::  F
end

SubzeroSet(f) = SubzeroSet{Float64}(f)
SubzeroSet{T}(f) where {T} = SubzeroSet{T,:closed}(f)
SubzeroSet{T,C}(f::F) where {T,C,F} = SubzeroSet{T,C,F}(f)

level(d::SubzeroSet) = 0

similardomain(d::SubzeroSet{S,C}, ::Type{T}) where {S,C,T} =
    SubzeroSet{T,C}(levelfun(d))

interior(d::SubzeroSet{T}) where {T} = SubzeroSet{T,:open}(d.f)
closure(d::SubzeroSet{T}) where {T} = SubzeroSet{T,:closed}(d.f)

boundary(d::SubzeroSet{T}) where {T} = ZeroSet{T}(d.f)

"The domain defined by `f(x) <= C` (or `f(x) < C`) for a given function `f` and constant `C`."
struct SublevelSet{T,C,F,S} <: AbstractSublevelSet{T,C}
    f       ::  F
    level   ::  S
end

SublevelSet(f, level) = SublevelSet{typeof(level)}(f, level)
SublevelSet{T}(f, level) where {T} = SublevelSet{T,:closed}(f, level)
SublevelSet{T,C}(f::F, level::S) where {T,C,F,S} = SublevelSet{T,C,F,S}(f, level)

similardomain(d::SublevelSet{S,C}, ::Type{T}) where {S,C,T} =
    SublevelSet{T,C}(levelfun(d), level(d))

interior(d::SublevelSet{T}) where {T} = SublevelSet{T,:open}(d.f, d.level)
closure(d::SublevelSet{T}) where {T} = SublevelSet{T,:closed}(d.f, d.level)

boundary(d::SublevelSet{T}) where {T} = LevelSet{T}(d.f, d.level)


"Supertype of superlevel set domains."
abstract type AbstractSuperlevelSet{T,C} <: FunctionLevelSet{T} end

indomain(x, d::AbstractSuperlevelSet{T,:closed}) where {T} = levelfun(d, x) >= level(d)
indomain(x, d::AbstractSuperlevelSet{T,:open}) where {T} = levelfun(d, x) > level(d)

show(io::IO, d::AbstractSuperlevelSet{T,:closed}) where {T} =
    print(io, "superlevel set f(x) >= $(level(d)) with f = $(levelfun(d))")
show(io::IO, d::AbstractSuperlevelSet{T,:open}) where {T} =
    print(io, "superlevel set f(x) > $(level(d)) with f = $(levelfun(d))")

==(d1::AbstractSuperlevelSet, d2::AbstractSuperlevelSet) =
    levelfun(d1)==levelfun(d2) && level(d1)==level(d2)
hash(d::AbstractSuperlevelSet, h::UInt) =
    hashrec("AbstractSuperlevelSet", levelfun(d), level(d), h)


"The domain where `f(x) >= 0` (or `f(x) > 0`)."
struct SuperzeroSet{T,C,F} <: AbstractSuperlevelSet{T,C}
    f   ::  F
end

SuperzeroSet(f) = SuperzeroSet{Float64}(f)
SuperzeroSet{T}(f) where {T} = SuperzeroSet{T,:closed}(f)
SuperzeroSet{T,C}(f::F) where {T,C,F} = SuperzeroSet{T,C,F}(f)

level(d::SuperzeroSet) = 0

similardomain(d::SuperzeroSet{S,C}, ::Type{T}) where {S,C,T} =
    SuperzeroSet{T,C}(levelfun(d))

interior(d::SuperzeroSet{T}) where {T} = SuperzeroSet{T,:open}(d.f)
closure(d::SuperzeroSet{T}) where {T} = SuperzeroSet{T,:closed}(d.f)

boundary(d::SuperzeroSet{T}) where {T} = ZeroSet{T}(d.f)


"The domain defined by `f(x) >= C` (or `f(x) > C`) for a given function `f` and constant `C`."
struct SuperlevelSet{T,C,F,S} <: AbstractSuperlevelSet{T,C}
    f       ::  F
    level   ::  S
end

SuperlevelSet(f, level) = SuperlevelSet{typeof(level)}(f, level)
SuperlevelSet{T}(f, level) where {T} = SuperlevelSet{T,:closed}(f, level)
SuperlevelSet{T,C}(f::F, level::S) where {T,C,F,S} = SuperlevelSet{T,C,F,S}(f, level)

similardomain(d::SuperlevelSet{S,C}, ::Type{T}) where {S,C,T} =
    SuperlevelSet{T,C}(levelfun(d), level(d))

interior(d::SuperlevelSet{T}) where {T} = SuperlevelSet{T,:open}(d.f, d.level)
closure(d::SuperlevelSet{T}) where {T} = SuperlevelSet{T,:closed}(d.f, d.level)

boundary(d::SuperlevelSet{T}) where {T} = LevelSet{T}(d.f, d.level)

## Additional functionality

pseudolevel(d::AbstractLevelSet, epsilon) =
    _pseudolevel(d, epsilon, levelfun(d), level(d))
_pseudolevel(d::AbstractLevelSet{T}, epsilon, fun, C) where {T} =
    SublevelSet{T,:open}(x -> norm(fun(x)-C), epsilon)

pseudolevel(d::Domain, epsilon) = pseudolevel(convert(LevelSet, d), epsilon)
