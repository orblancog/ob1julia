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
        for k in 1:4
            @inbounds pp = pp+A[i][k]*B[k][j]
        end
        push!(p,pp)
    end
    q::Array{Basic,2}=transpose(reshape(p,4,4))
    return q
end
function dotransport(mlist)
    for i in 1:length(mlist)-1
        if i == 1
            global m1 = mlist[1]
        else
            global m1 = mult(m1,mlist[i]) ### mult modifies mt
        end
    end
    return m1
end

### BEGIN of the script

println("  OB1's APOCHROMAT DESIGN 0.1")
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

m = dotransport(mlist)
println(typeof(m[1][1]))
println(expand(m[1][1]))
dm_dE=diff(expand(m[1][1]),dE)
d2m_dE2=diff(dm_dE,dE)
f0 = subs(expand(m[1][1]),dE,0)
print(f0)
#mytaylor(expand[1][1])
exit()


exit()
println("sub0        : ",subs(m[1][1],dE => 0))


        
println("dm_dE * dE  :",dm_dE*dE)
println("d2m_dE2 * dE:",d2m_dE2*dE^2)



exit()


dRq_12_dE = diff(Rq_21,dE)
println(dRq_12_dE)
exit()

Rd = [symbols("Rd_$i$j") for i in 1:4, j in 1:4]
Rd11 = 1
Rd12 = d
Rd21 = 0
Rd22 = 1

R  = [symbols("R_$i$j") for i in 1:4, j in 1:4]
println(R)
for i in 1:4, j in 1:4
    R_ij = R_ij + Rqij*Rqdij
end
exit()
println(expand(Rm[1][1]))
exit()
println(diff(Rm[1][1],dE))
#println(diff(Rm[1][2],dE))

#println(diff(Rm[2][1],dE))
#println(diff(Rm[2][2],dE))
#println(length(Rm))
