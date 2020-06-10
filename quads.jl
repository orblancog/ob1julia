### This script is used to calculate the quads gradients k, lengths l,
###     and drifts d for an apochromatic transport line
### author : orblancog@gmail.com
### 2020/06/07

### Written in Julia 1.5

#### usign algebras
using SymEngine
using LinearAlgebra

### energy parameter
dE=symbols(:dE);

### structs to create the elements
struct quad
    k :: Basic
    l :: Basic
end
struct drift
    l :: Basic
end
### Matrix declaration
function mq(q)
    # horizontally  focusing quadrupole transport matrix
    # four rows
    # quad chromaticity k->k(1+dE)
    R=[Basic[1 0 0 0];Basic[-1/(q.k*q.l*(1+dE)) 1 0 0];Basic[0 0 1 0];Basic[0 0 1/(q.k*(1+dE)*q.l) 1]]
    return R
end
function md(d)
    # four rows
    R=[Basic[1 d.l 0 0];Basic[0 1 0 0];Basic[0 0 1 d.l];Basic[0 0 0 1]]
    return R
end

### Matrix multiplication
function mult(A,B)
    p = Any[]
    for i in 1:4, j in 1:4
        pp = 0
        typeof(pp)
        for k in 1:4
            pp = pp + A[i,k]*B[k,j]### order changed because push is adding column vectors
#            println("i",i," j",j," k",k," pp",pp," A",A[i,k]," B",B[k,j])
        end
        push!(p,pp) ### maybe push is adding column vectors
#        println(" stato p",p)
    end
    return transpose(reshape(p,4,4))
end
function dotransport(mlist)
    global m1 = mlist[1]
    for i in 1:length(mlist)-1
        if length(mlist) > 1
            #            println("enter")
            global m1 = mult(m1,mlist[i+1]) ### mult modifies mt
        end
    end
    return m1
end
### second order taylor expansion
function dotaylor(p)
#    print("create")
    f   = expand(p)
    df  = diff(f,  dE)
    ddf = diff(df, dE)
    ff = [symbols("f$i") for i in 0:1:2]
    f0 = subs(f  , dE, 0)
    f1 = subs(df , dE, 0)
    f2 = subs(ddf, dE, 0)

    ftaylor = expand(f0 + f1*dE + 0.5*f2*dE^2)
    return ftaylor
end
### truncate expansion
function dotruncate(p)
    pexpnd = expand(p)
    ptrunc = 0
    for i in 0:1:2 ## zeroth to second order expansion
        ptrunc = ptrunc + coeff(pexpnd,dE,Basic(i))*dE^i
    end
    return expand(ptrunc)
end
### BEGIN of the script

println("  OB1's APOCHROMAT DESIGN 0.2")
Nq = 2 # number of quadrupoles
Nl = Nq - 1 # number of drifts
println("    Using ",Nq," quads and ",Nl," drifts")

### Create variables
println("    Creating symbols...")
#R = [symbols("R$i$j") for i in 1:4, j in 1:4]
kk = [symbols("k$i") for i in 1:Nq]
ll = [symbols("l$i") for i in 1:Nq]
dd = [symbols("d$i") for i in 1:Nl]
lstar = symbols("lstar")
println("    symbols created.")

## create elements
println("    Creating elements...")
dstar = drift(lstar)
q1 = quad(kk[1],ll[1])
d1 = drift(dd[1])
println("    elements created.")

### create matrices
println("    Creating matrix representation of the transport line...")
mlist = Any[]
mdstar = md(dstar)
push!(mlist,mdstar) # last element first
push!(mlist, mq(quad(kk[Nq],ll[Nq])))
for i in Nl:-1:1
    push!(mlist, md(drift(dd[i])))
    push!(mlist, mq(quad(kk[i],ll[i])))
end
push!(mlist,mdstar) # first elemen last 
println("    ",length(mlist)," matrices created.")
println("    Transport line created. Adios !")

println("    Creating the transport matrix to second order in dE")
#m = dotransport(mlist)

mlist2 = Any[]
#push!(mlist2,mdstar)
push!(mlist2,mdstar)
push!(mlist2,mq(quad(kk[1],ll[1])))
#push!(mlist2,mq(quad(kk[2],ll[2])))
#push!(mlist2,mdstar)
#println("length",length(mlist2))
R = dotransport(mlist2)
println("    matrix ",R)
Rarraytaylor = Any[]
for i in 1:4, j in 1:4
    Rij = dotaylor(expand(R[i,j]))
    push!(Rarraytaylor,Rij)
end
Rtaylor=transpose(reshape(Rarraytaylor,4,4))
println("    Taylor expansion to second order of the matrix R")
for i in 1:4, j in 1:4
    println("      R$i$j = ",expand(Rtaylor[i,j]))
end

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
println("      betaxtrunc = ",betaxtrunc)
println("      betaytrunc = ",betaytrunc)
println("      alfaxtrunc = ",alfaxtrunc)
println("      alfaytrunc = ",alfaytrunc)
println("      gamaxtrunc = ",gamaxtrunc)
println("      gamaytrunc = ",gamaytrunc)
exit()
