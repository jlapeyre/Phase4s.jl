module Phase4s

using Random: Random

export Phase, phase_from_factors

# UInt8: minimum storage req. Best performance in all of a few tests.
const PT = UInt8

"""
    struct Phase

Implements the cyclic group ℤ₄ represented by the fourth roots of unity.

Conversion and compatibility with other number types is supported.

# Examples:

Construction from numbers, such as `Int`, `Complex{Int}`, `Float64`, etc:
```jl-doctest
julia> Phase.((1, im, -1, -im))
(Phase(+1), Phase(+im), Phase(-1), Phase(-im))

julia> Phase.((1., im * 1.0, -1.))
(Phase(+1), Phase(+im), Phase(-1))

julia> Phase(2)
ERROR: ArgumentError: Can't convert `2` to a Pauli phase
```

Construction with values of the underlying storage type, `UInt8`, is treated differently:
```jl-doctest
julia> Phase.((0x00, 0x01, 0x02, 0x03))
(Phase(+1), Phase(+im), Phase(-1), Phase(-im))

julia> Phase(0x04)
ERROR: ArgumentError: `0x04` not a valid Pauli phase.
```

Multiplication, power, inverse:
```jl-doctest
julia> Phase(im) .* Phase.((1, im, -1, -im))
(Phase(+im), Phase(-1), Phase(-im), Phase(+1))

julia> inv.(Phase.((1, im, -1, -im)))
(Phase(+1), Phase(-im), Phase(-1), Phase(+im))

julia> Phase(im)^5
Phase(+im)
```

Multiplication by other number types:
```jl-doctest
julia> (1 * Phase(im), 1.0 * Phase(im), Phase(-im) * 2.1)
(0 + 1im, 0.0 + 1.0im, 0.0 - 2.1im)
```

Other functions:
```jl-doctest
julia> (one(Phase), one(Phase(im)))
(Phase(+1), Phase(+1))

julia> x = convert(Complex{Int128}, Phase(im)); (x, typeof(x))
(0 + 1im, Complex{Int128})
```
"""
struct Phase <: Number
    x::PT
    function Phase(x::PT)
        x < PT(4) ||
            throw(ArgumentError(lazy"`$(_string_uint(x))` not a valid Pauli phase."))
        return new(x)
    end

    function Phase(x::Union{Real,Complex})
        x == one(x) && return Phase(PT(0))
        x == Complex(0, 1) && return Phase(PT(1))
        x == -one(x) && return Phase(PT(2))
        x == -Complex(0, 1) && return Phase(PT(3))
        throw(ArgumentError(lazy"Can't convert `$x` to a Pauli phase"))
    end
end

function _string_uint(x::Unsigned)
    return "0x" * string(x; pad=sizeof(x) << 1, base=16)
end
# Use `im` rather than `i` so that we print a constructor call.
function _print_phase(io::IO, p::Phase)
    return print(io, typeof(p), "(", ("+1", "+im", "-1", "-im")[p.x + 1], ")")
end
Base.show(io::IO, ::MIME"text/plain", p::Phase) = _print_phase(io, p)
Base.show(io::IO, p::Phase) = _print_phase(io, p)

Base.one(::Type{Phase}) = Phase(1)
Base.one(::Phase) = one(Phase)
# Two reasons for these MethodError's
# 1) To add a method, someone will have to actively remove these methods, and hopefully
#    ask why they are here.
# 2) The location of the throw, and the error message is a bit more understandable compared
#    to throwing in the fallback methods.
Base.iszero(p::Phase) = throw(MethodError(iszero, (p,)))
Base.iszero(::Type{Phase}) = throw(MethodError(iszero, (Phase,)))
Base.zero(p::Phase) = throw(MethodError(zero, (p,)))
Base.zero(::Type{Phase}) = throw(MethodError(zero, (Phase,)))

Base.complex(p::Phase) = p
Base.float(p::Phase) = Complex{Float64}(p)

# This is as fast as checking the low bit explicitly v1.11.2
Base.isreal(p::Phase) = p == Phase(1) || p == Phase(-1)

Base.isinteger(p::Phase) = isreal(p)
Base.real(p::Phase) = real(Int, p)
Base.real(::Type{T}, p::Phase) where {T} = reim(T, p)[1]
Base.imag(p::Phase) = imag(Int, p)
Base.imag(::Type{T}, p::Phase) where {T} = reim(T, p)[2]

# Faster than if-else v1.11.2
function Base.reim(::Type{T}, p::Phase) where {T}
    vals = ((one(T), zero(T)), (zero(T), one(T)), (-one(T), zero(T)), (zero(T), -one(T)))
    return @inbounds vals[p.x + 1]
end
Base.reim(p::Phase) = reim(Int, p)
Base.ispow2(p::Phase) = p == Phase(1)
# Faster than if-else v1.11.2
Base.sign(p::Phase) = @inbounds (Phase(1), Phase(1), Phase(-1), Phase(-1))[p.x + 1]

# Don't need this method. It already throws immediately
# Do not define a method for `signbit`; Base does not define it for Complex.
# Base.signbit(p::Phase) = throw(MethodError(signbit, (p,)))

# Fallback would be ok, if Phase <: AbstractComplex. But AbstractComplex does not exist.
Base.flipsign(p::Phase, y::Real) = ifelse(signbit(y), -p, p)
Base.in(p::Phase, r::AbstractRange{<:Real}) = isreal(p) && real(p) in r

Base.:(*)(l::Phase, r::Phase) = Phase(mod(l.x + r.x, PT(4)))
Base.:(/)(l::Phase, r::Phase) = Phase(mod(l.x - r.x, PT(4)))

Base.convert(T::Type{<:Complex}, p::Phase) = (one(T), T(im), T(-one(T)), T(-im))[p.x + 1]
# Ensure a concrete type for Complex(::Phase)
Base.Complex(p::Phase) = Complex{Int}(p)
(::Type{T})(p::Phase) where {T<:Complex} = convert(T, p)
Base.promote_rule(::Type{Phase}, T::Type{<:Complex}) = T
Base.promote_rule(::Type{Phase}, ::Type{Complex{Bool}}) = promote_rule(Phase, Complex{Int})
Base.promote_rule(::Type{Phase}, T::Type{<:Real}) = Complex{T}
Base.promote_rule(::Type{Phase}, ::Type{Bool}) = promote_rule(Phase, Complex{Int})
# For *(::Phase, ::Number) and /, the fallbacks are efficient.

# Unary minus
Base.:(-)(p::Phase) = Phase(mod(p.x + PT(2), PT(4)))
# Power
Base.:(^)(p::Phase, n::Integer) = Phase(mod(n * p.x, PT(4)) % PT)
Base.:(^)(p::Phase, b::Bool) = ifelse(b, p, Phase(1))
Base.inv(p::Phase) = Phase(mod(-p.x, PT(4)))

# `Phase` is not closed under the following operations
for op in (:(+), :(-))
    @eval Base.$op(p1::Phase, p2::Phase) = $op(Complex(p1), Complex(p2))
end

for op in (:(^),)
   # No red in @code_warntype on v1.11.2
   @eval let  v = [complex(float(x)) for x in (1,im,-1,-im)], res = [$op(x,y) for x in v, y in v]
        Base.$op(p1::Phase, p2::Phase) = @inbounds res[p1.x + 4 * p2.x + 1]
    end
end

Base.conj(p::Phase) = iseven(p.x) ? p : -p
Base.log(p::Phase) = Complex{Float64}(0, angle(p))
Base.abs(p::Phase) = 1
Base.abs2(p::Phase) = 1

for op in (:exp, :angle, :cos, :sin, :tan, :asin, :acos, :atan,
           :cosh, :sinh, :tanh, :acosh, :asinh, :atanh,
           :cis, :cispi, :log2, :log10, :expm1, :log1p, :exp10, :exp2)
   @eval Base.$op(p::Phase) = ($op(1.0 + 0im), $op(0.0 + 1.0im), $op(-1.0 + 0im), $op(0.0 - 1.0im))[p.x + 1]
end

"""
    phase_from_factors(num_imag, num_minus)::Phase

Construct a Pauli phase from `num_imag` factors of the imaginary unit
and `num_minus` factors of `-1`

This is equivalent to `Phase(im^num_imag * (-1)^num_minus)`.
"""
function phase_from_factors(num_imag, num_minus)::Phase
    # To multiply by -1, rotate two elements to the right
    return Phase(PT(mod(num_imag + (!iseven(num_minus) << 1), PT(4))))
end

# A bit faster than rand(0x00:0x03) v1.11.2
function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Phase})
    # Following has optimization possibility, but is actually slower in 1.11.
    # (a, b) = rand(Tuple{Bool,Bool})
    # Following is slightly faster v1.11.2
    (a, b) = (rand(Bool), rand(Bool))

    # Converting a and b to UInt8 separately is faster v1.11.2
    return Phase(PT(a) | PT(b) << 1)
end

end # module Phase4s
