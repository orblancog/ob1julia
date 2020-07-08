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
println("  OB1's APOCHROMAT DESIGN 2.0")
println("     author  : orblancog@gmail.com ")
println("     release : 2020/jun/30")
println("")
println("    Creating the transport matrix to second order in dE")
#R = dotransport(mlist2)
R = dotransport(mpt104())
#R = dotransport(mlist4)
#copympt = mpt104v2D()
#copympt[end] = [1 0; 0 1]
#println(copympt[end])
#R = dotransport(copympt)

println("    Taylor expansion to second order of the matrix R")
### taylor expansion already done during multiplication
Rtaylor = R

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
Sx  = Rtaylor[1,2]
Cpx = Rtaylor[2,1]
Spx = Rtaylor[2,2]
Cy  = Rtaylor[3,3]
Sy  = Rtaylor[3,4]
Cpy = Rtaylor[4,3]
Spy = Rtaylor[4,4]

betax =   (Cx^2)*betax0  -           2*Cx*Sx*alfax0 +  (Sx^2)*gamax0
alfax =  -Cx*Cpx*betax0  + (Cx*Spx + Sx*Cpx)*alfax0 -  Sx*Spx*gamax0
gamax =  (Cpx^2)*betax0  -         2*Cpx*Spx*alfax0 + (Spx^2)*gamax0
betay =   (Cy^2)*betay0  -           2*Cy*Sy*alfay0 +  (Sy^2)*gamay0
alfay =  -Cy*Cpy*betay0  + (Cy*Spy + Sy*Cpy)*alfay0 -  Sy*Spy*gamay0
gamay =  (Cpy^2)*betay0  -         2*Cpy*Spy*alfay0 + (Spy^2)*gamay0

println("    ... twiss parameters defined")

println("    Expanding twiss parameters to second order in dE")
betaxtrunc = dotaylor(betax)
alfaxtrunc = dotaylor(alfax)
gamaxtrunc = dotaylor(gamax)
betaytrunc = dotaylor(betay)
alfaytrunc = dotaylor(alfay)
gamaytrunc = dotaylor(gamay)


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
#=
eqbetaxdE0=betaxtrunc
eqbetaxdE1=diff(betaxtrunc,dE)*dE
eqbetaxdE2=0.5*diff(diff(betaxtrunc,dE),dE)*dE^2
eqbetaydE0=betaytrunc
eqbetaydE1=diff(betaytrunc,dE)*dE
eqbetaydE2=0.5*diff(diff(betaytrunc,dE),dE)*dE^2
eqalfaxdE0=alfaxtrunc
eqalfaxdE1=diff(alfaxtrunc,dE)*dE
eqalfaxdE2=0.5*diff(diff(alfaxtrunc,dE),dE)*dE^2
eqalfaydE0=alfaytrunc
eqalfaydE1=diff(alfaytrunc,dE)*dE
eqalfaydE2=0.5*diff(diff(alfaytrunc,dE),dE)*dE^2
=#

eqbetaxdE0=subs(betaxtrunc, dE, 0)
eqbetaxdE1=subs(diff(betaxtrunc,dE), dE, 0 )
eqbetaxdE2=0.5*subs(diff(diff(betaxtrunc,dE),dE), dE, 0)
eqalfaxdE0=subs(alfaxtrunc, dE, 0)
eqalfaxdE1=subs(diff(alfaxtrunc,dE), dE, 0)
eqalfaxdE2=0.5*subs(diff(diff(alfaxtrunc,dE),dE), dE, 0)


eqbetaydE0=subs(betaytrunc, dE, 0)
eqbetaydE1=subs(diff(betaytrunc,dE), dE, 0)
eqbetaydE2=0.5*subs(diff(diff(betaytrunc,dE),dE), dE, 0)
eqalfaydE0=subs(alfaytrunc, dE, 0)
eqalfaydE1=subs(diff(alfaytrunc,dE), dE, 0)
eqalfaydE2=0.5*subs(diff(diff(alfaytrunc,dE),dE), dE, 0)

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


#x
exprbetx0=expand(subs(subs(subs(eqbetaxdE0, alfax0, 0), gamax0, 1/betax0),betax0, 0.2))
exprbetx1=expand(subs(subs(subs(eqbetaxdE1, alfax0, 0), gamax0, 1/betax0),betax0, 0.2))
exprbetx2=expand(subs(subs(subs(eqbetaxdE2, alfax0, 0), gamax0, 1/betax0),betax0, 0.2))
expralfx0=expand(subs(subs(subs(eqalfaxdE0, alfax0, 0), gamax0, 1/betax0),betax0, 0.2))
expralfx1=expand(subs(subs(subs(eqalfaxdE1, alfax0, 0), gamax0, 1/betax0),betax0, 0.2))
expralfx2=expand(subs(subs(subs(eqalfaxdE2, alfax0, 0), gamax0, 1/betax0),betax0, 0.2))
#println("      eqbetaxdE0 = ",exprbetx0)
#println("      eqbetaxdE1 = ",exprbetx1)
#println("      eqbetaxdE2 = ",exprbetx2)
#println("      eqalfaxdE0 = ",expralfx0)
#println("      eqalfaxdE1 = ",expralfx1)
#println("      eqalfaxdE2 = ",expralfx2)
apolog=open("apo.log","w")
#write(apolog,"betax(dE) = ",exprbetx0," + 1*",exprbetx1,"*dE + 1*",exprbetx2,"*dE**2")
#write(apolog,"alfax(dE) = ",expralfx0," + 1*",expralfx1,"*dE + 1*",expralfx2,"*dE**2")
write(apolog,"betax(dE) = $exprbetx0 + 1*$exprbetx1*dE + 1*$exprbetx2*dE**2\n")
write(apolog,"alfax(dE) = $expralfx0 + 1*$expralfx1*dE + 1*$expralfx2*dE**2\n")
write(apolog,"")


#y
exprbety0=expand(subs(subs(subs(eqbetaydE0, alfay0, 0), gamay0, 1/betay0),betay0, 0.2))
exprbety1=expand(subs(subs(subs(eqbetaydE1, alfay0, 0), gamay0, 1/betay0),betay0, 0.2))
exprbety2=expand(subs(subs(subs(eqbetaydE2, alfay0, 0), gamay0, 1/betay0),betay0, 0.2))
expralfy0=expand(subs(subs(subs(eqalfaydE0, alfay0, 0), gamay0, 1/betay0),betay0, 0.2))
expralfy1=expand(subs(subs(subs(eqalfaydE1, alfay0, 0), gamay0, 1/betay0),betay0, 0.2))
expralfy2=expand(subs(subs(subs(eqalfaydE2, alfay0, 0), gamay0, 1/betay0),betay0, 0.2))
#=
println("      eqbetaydE0 = ",exprbety0)
println("      eqbetaydE1 = ",exprbety1)
println("      eqbetaydE2 = ",exprbety2)
println("      eqalfaydE0 = ",expralfy0)
println("      eqalfaydE1 = ",expralfy1)
println("      eqalfaydE2 = ",expralfy2)
=#
write(apolog,"betay(dE) = $exprbety0 + 1*$exprbety1*dE + 1*$exprbety2*dE**2\n")
write(apolog,"alfay(dE) = $expralfy0 + 1*$expralfy1*dE + 1*$expralfy2*dE**2\n")
close(apolog)
println(" file  apo.log created")
println(" Chao, 8D !!!")
exit()

