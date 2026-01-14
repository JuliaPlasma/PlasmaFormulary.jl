# Workaround for UnitfulEquivalences.uconvert not supporting conversions of
# quantites with the same dimensions. See
# https://github.com/sostock/UnitfulEquivalences.jl/issues/19
function custom_uconvert(dest_unit, x, equivalence)
    return if dimension(x) == dimension(dest_unit)
        # plain Unitful ustrip
        Unitful.uconvert(dest_unit, x)
    else
        UnitfulEquivalences.uconvert(dest_unit, x, equivalence)
    end
end

_to_Hz(ω) = uconvert(u"Hz", ω, Periodic())
_to_Hz(ω, flag) = flag ? _to_Hz(ω) : ω

function _norm(𝐱)
    @assert length(𝐱) == 3
    return sqrt(norm_sqr(𝐱[1]) + norm_sqr(𝐱[2]) + norm_sqr(𝐱[3]))
end
