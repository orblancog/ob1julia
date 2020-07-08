### This script is used to calculate the quads gradients k, lengths l,
###     and drifts d for an apochromatic transport line
### author : orblancog@gmail.com
### 2020/06/07

### Written in Julia 1.5

#### usign algebras
using SymEngine
using LinearAlgebra

## first order apochromat, mpt104
include("mpt104.jl")
include("mpt104v2D.jl")
include("elements.jl")
include("myfunctions.jl")

### energy parameter
dE=symbols(:dE)

### BEGIN of the script
println("  OB1's APOCHROMAT DESIGN 2.2")
println("     author  : orblancog@gmail.com ")
println("     release : 2020/jun/30")
println("")
println("    Creating the transport matrix to second order in dE")
#R = dotransport(mlist2)

#R = dotransport(mlist4)
#copympt = mpt104v2D()
#copympt[end] = [1 0; 0 1]
#println(copympt[end])
#R = dotransport(copympt)

for i in 1:2
    @show i
    local R = dotransport(mpt104())
    twisspropagation(R,0.2,0.2,0,0,1/0.2,1/0.2)
    
end

println(" Chao, 8D !!!")
exit()

