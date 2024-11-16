using Test
using PlasmaFormulary: @permutable_args

# Test function to verify permutable arguments
@permutable_args function test_func(x::Int, y::String)
    return "x=$x, y=$y"
end

@testset "permutable_args macro" begin
    # Test normal order
    @test test_func(42, "hello") == "x=42, y=hello"

    # Test reversed order
    @test test_func("hello", 42) == "x=42, y=hello"

    # Test type errors
    @test_throws MethodError test_func(1.0, "hello")  # Wrong type for x
    @test_throws MethodError test_func("hello", "world")  # Wrong type for y

    # Test that both methods exist
    methods_count = length(methods(test_func))
    @test methods_count == 2  # Should have generated 2 methods
end

@permutable_args function test_func_kw(x::Int, y::String; z::Float64 = 3.14)
    return "x=$x, y=$y, z=$z"
end

@testset "permutable_args macro with keyword arguments" begin
    @test test_func_kw(42, "hello") == "x=42, y=hello, z=3.14"
    @test test_func_kw(42, "hello", z = 2.71) == "x=42, y=hello, z=2.71"
    @test test_func_kw("hello", 42; z = 520.0) == "x=42, y=hello, z=520.0"
end

@testset "permutable_args macro - one-line function definition" begin
    @permutable_args test_func_one_line(x::Int, y::String; z::Float64 = 3.14) =
        "x=$x, y=$y, z=$z"
    @test test_func_one_line(42, "hello") == "x=42, y=hello, z=3.14"
    @test test_func_one_line(42, "hello", z = 2.71) == "x=42, y=hello, z=2.71"
    @test test_func_one_line("hello", 42; z = 520.0) == "x=42, y=hello, z=520.0"
end