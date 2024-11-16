using BenchmarkTools

test_func_base(x::Int, y::String, z::Float64, s::Symbol, b::Bool) =
    "x=$x, y=$y, z=$z, s=$s, b=$b"

@permutable_args function test_func_multiple(
    x::Int,
    y::String,
    z::Float64,
    s::Symbol,
    b::Bool,
)
    return test_func_base(x, y, z, s, b)
end

@btime test_func_base(1, "hello", 3.14, :s, false)
@btime test_func_multiple(1, "hello", 3.14, :s, false)
@btime test_func_multiple(false, "hello", 1, :s, 3.14)
