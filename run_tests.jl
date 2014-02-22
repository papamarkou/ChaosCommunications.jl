examples = ["UTURUniChannelCSK_example"]

println("Running tests:")

for t in examples
    test_fn = joinpath("examples", "$t.jl")
    println("  * $test_fn *")
    include(test_fn)
end
