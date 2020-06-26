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
include("elements.jl")
include("myfunctions.jl")

### energy parameter
dE=symbols(:dE);


### BEGIN of the script

println("  OB1's APOCHROMAT DESIGN 0.6")
Nq = 10 # number of quadrupoles
Nl = Nq - 1 # number of drifts
#println("    Using ",Nq," quads and ",Nl," drifts")

### Create variables
println("    Creating symbols...")
#R = [symbols("R$i$j") for i in 1:4, j in 1:4]
kk = [symbols("k$i") for i in 1:Nq]
ll = [symbols("l$i") for i in 1:Nq]
kl = [symbols("kl$i") for i in 1:Nq]
dd = [symbols("d$i") for i in 1:Nl]
lstar = symbols("lstar")
println("    symbols created.")

## create elements
println("    Creating elements...")
dstar = drift(lstar)
#q1  = quad(kk[1],ll[1])
#qq1 = qquad(kl[1])#other type of quad, I hope faster
#d1  = drift(dd[1])
println("    elements created.")

### create matrices
println("    Creating matrix representation of the transport line...")
mlist = Any[]
#mdstar = md(dstar)
mdstar = md(lstar)
push!(mlist,mdstar) # last element first
#push!(mlist, mq(quad(kk[Nq],ll[Nq])))
push!(mlist, mqq(qquad(kl[Nq])))
for i in Nl:-1:1
#    push!(mlist, md(drift(dd[i])))
    push!(mlist, md(1.0))
#    push!(mlist, mq(quad(kk[i],ll[i])))
    push!(mlist, mqq(qquad(kl[i])))
end
push!(mlist,mdstar) # first elemen last 
println("    ",length(mlist)," matrices created.")
println("    Transport line created. Adios !")

println("    Creating the transport matrix to second order in dE")
#m = dotransport(mlist)

mlist2 = Any[]
push!(mlist2,mq(quad(kk[1],ll[1])))
push!(mlist2,mdstar)
#push!(mlist2,mdstar)
# #push!(mlist2,mq(quad(kk[2],ll[2])))
# #push!(mlist2,mdstar)
# println("length",length(mlist2))
# println(mlist2)
# exit()
#R = dotransport(mlist2)

R = dotransport(mpt104())

#=
#println("    matrix ",R)
for i in 1:4, j in 1:4
    println("      R$i$j = ",R[i,j])
end
=#



Rarraytaylor = Any[]
for i in 1:4, j in 1:4
#        if i==1 && j==2
            Rij = dotaylor(expand(R[i,j]))
            push!(Rarraytaylor,Rij)
#        end
end
#exit()

Rtaylor=transpose(reshape(Rarraytaylor,4,4))

println("    Taylor expansion to second order of the matrix R")
#=
for i in 1:4, j in 1:4

        println("      R$i$j = ",expand(Rtaylor[i,j]))
    end
end
=#
#exit()
println("    Calculation of twiss parameters...")
# twiss symbols as a function of dE
betax = symbols("betax")
betay = symbols("betay")
alfax = symbols("alfax")
alfay = symbols("alfay")
gamax = symbols("gamax")
gamay = symbols("gamay")
# initial values
betax0 = symbols("betax0")
betay0 = symbols("betay0")
alfax0 = symbols("alfax0")
alfay0 = symbols("alfay0")
gamax0 = symbols("gamax0")
gamay0 = symbols("gamay0")
# notation to propagate twiss
Cx = symbols("Cx")
Cy = symbols("Cy")
Sx = symbols("Sx")
Sy = symbols("Sy")
Cx  = Rtaylor[1,1]
Cy  = Rtaylor[3,3]
Sx  = Rtaylor[1,2]
Sy  = Rtaylor[3,4]
Cpx = Rtaylor[2,1]
Cpy = Rtaylor[4,3]
Spx = Rtaylor[2,2]
Spy = Rtaylor[4,4]

betax =  (Cx^2)*betax0  - 2*Cx*Sx*alfax0 + (Sx^2)*gamax0
betay =  (Cy^2)*betay0  - 2*Cy*Sy*alfay0 + (Sy^2)*gamay0
alfax = -Cx*Cpx*betax0 + (Cx*Spx + Sx*Cpx)*alfax0 - Sx*Spx*gamax0
alfay = -Cy*Cpy*betay0 + (Cy*Spy + Sy*Cpy)*alfay0 - Sy*Spy*gamay0
gamax =  (Cpx^2)*betax0 - 2*Cpx*Spx*alfax0 + (Spx^2)*gamax0
gamay =  (Cpy^2)*betay0 - 2*Cpy*Spy*alfay0 + (Spy^2)*gamay0
println("    ... twiss parameters defined")

println("    Expanding twiss parameters to second order in dE")
betaxtrunc = dotruncate(betax)
betaytrunc = dotruncate(betay)
alfaxtrunc = dotruncate(alfax)
alfaytrunc = dotruncate(alfay)
gamaxtrunc = dotruncate(gamax)
gamaytrunc = dotruncate(gamay)
#=
println("      betaxtrunc = ",betaxtrunc)
println("      betaytrunc = ",betaytrunc)
println("      alfaxtrunc = ",alfaxtrunc)
println("      alfaytrunc = ",alfaytrunc)
println("      gamaxtrunc = ",gamaxtrunc)
println("      gamaytrunc = ",gamaytrunc)
=#
println("    Expansion finished")
#exit()
println("    Extracting chromatics derivatives")
eqbetaxdE0=betaxtrunc
eqbetaxdE1=diff(betaxtrunc,dE)*dE
eqbetaxdE2=0.5diff(diff(betaxtrunc,dE),dE)*dE^2
eqbetaydE0=betaytrunc
eqbetaydE1=diff(betaytrunc,dE)*dE
eqbetaydE2=0.5diff(diff(betaytrunc,dE),dE)*dE^2
eqalfaxdE0=alfaxtrunc
eqalfaxdE1=diff(alfaxtrunc,dE)*dE
eqalfaxdE2=0.5diff(diff(alfaxtrunc,dE),dE)*dE^2
eqalfaydE0=alfaytrunc
eqalfaydE1=diff(alfaytrunc,dE)*dE
eqalfaydE2=0.5diff(diff(alfaytrunc,dE),dE)*dE^2
#=
println("      eqbetaxdE0 = ",eqbetaxdE0)
println("      eqbetaxdE1 = ",eqbetaxdE1)
println("      eqbetaxdE2 = ",eqbetaxdE2)
println("      eqbetaydE0 = ",eqbetaydE0)
println("      eqbetaydE1 = ",eqbetaydE1)
println("      eqbetaydE2 = ",eqbetaydE2)
println("      eqalfaxdE0 = ",eqalfaxdE0)
println("      eqalfaxdE1 = ",eqalfaxdE1)
println("      eqalfaxdE2 = ",eqalfaxdE2)
println("      eqalfaydE0 = ",eqalfaydE0)
println("      eqalfaydE1 = ",eqalfaydE1)
println("      eqalfaydE2 = ",eqalfaydE2)
=#
println("      eqbetaxdE0 = ",expand(subs(subs(eqbetaxdE0, alfax0, 0), gamax0, 1/betax0)))
#println("      eqbetaxdE1 = ",expand(subs(subs(eqbetaxdE1, alfax0, 0), gamax0, 1/betax0)))
#println("      eqbetaxdE2 = ",expand(subs(subs(eqbetaxdE2, alfax0, 0), gamax0, 1/betax0)))
println("      eqbetaydE0 = ",expand(subs(subs(eqbetaydE0, alfay0, 0), gamay0, 1/betay0)))
#println("      eqbetaydE1 = ",expand(subs(subs(eqbetaydE1, alfay0, 0), gamax0, 1/betay0)))
#println("      eqbetaydE2 = ",expand(subs(subs(eqbetaydE2, alfay0, 0), gamax0, 1/betay0)))
#println("      eqalfaxdE0 = ",expand(subs(subs(eqalfaxdE0, alfax0, 0), gamax0, 1/betax0)))
#println("      eqalfaxdE1 = ",expand(subs(subs(eqalfaxdE1, alfax0, 0), gamax0, 1/betax0)))
#println("      eqalfaxdE2 = ",expand(subs(subs(eqalfaxdE2, alfax0, 0), gamax0, 1/betax0)))
#println("      eqalfaydE0 = ",expand(subs(subs(eqalfaydE0, alfax0, 0), gamax0, 1/betax0)))
#println("      eqalfaydE1 = ",expand(subs(subs(eqalfaydE1, alfax0, 0), gamax0, 1/betax0)))
#println("      eqalfaydE2 = ",expand(subs(subs(eqalfaydE2, alfax0, 0), gamax0, 1/betax0)))


exit()
