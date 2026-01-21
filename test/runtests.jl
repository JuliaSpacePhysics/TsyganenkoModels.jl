using TsyganenkoModels
using Test
using Aqua

@testset "TsyganenkoModels.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(TsyganenkoModels)
    end
    # Write your tests here.
end
