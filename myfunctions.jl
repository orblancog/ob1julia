
### Matrix declaration
### deprecated thin quad
#=
function mq(q)
    # horizontally  focusing quadrupole transport matrix
    # four rows
    # quad chromaticity k->k/(1+dE)
    R=[Basic[1 0 0 0]; Basic[-(1+dE)/(q.k*q.l) 1 0 0]; Basic[0 0 1 0];Basic[0 0 (1+dE)/(q.k*q.l) 1]]
    return R
end
=#
# focusing quad
function mqf(q)
    # horizontally  focusing quadrupole transport matrix
    # four rows
    # quad chromaticity k->k/(1+dE)
    Rf11 =   cos(q.l*sqrt(abs(q.k)/(1+dE)))
    Rf12 =  sin(q.l*sqrt(abs(q.k)/(1+dE))) / (sqrt(abs(q.k)/(1+dE)))
    Rf21 = -(sqrt(abs(q.k)/(1+dE)))*sin(q.l*sqrt(abs(q.k)/(1+dE)))
    Rf22 =  Rf11
    println("      quad foc")
    R=[Basic[Rf11 Rf12]; Basic[Rf21 Rf22]]
    return R
end
# defocusing quad
function mqd(q)
    # horizontally  focusing quadrupole transport matrix
    # four rows
    # quad chromaticity k->k/(1+dE)
    Rd11 =  cosh(q.l*sqrt(abs(q.k)/(1+dE)))
    Rd12 =  sinh(q.l*sqrt(abs(q.k)/(1+dE)))/(sqrt(abs(q.k)/(1+dE)))
    Rd21 =  (sqrt(abs(q.k)/(1+dE)))*sinh(q.l*sqrt(abs(q.k)/(1+dE)))
    Rd22 =  Rd11
    println("      quad defoc")
    R=[Basic[Rd11 Rd12]; Basic[Rd21 Rd22]]
    return R
end

# thick quad 2D
function mq(q)
    # horizontally  focusing quadrupole transport matrix
    # four rows
    # quad chromaticity k->k/(1+dE)
    Rf11 =   cos(q.l*sqrt(abs(q.k)/(1+dE)))
    Rf12 =  sin(q.l*sqrt(abs(q.k)/(1+dE))) / (sqrt(abs(q.k)/(1+dE)))
    Rf21 = -(sqrt(abs(q.k)/(1+dE)))*sin(q.l*sqrt(abs(q.k)/(1+dE)))
    Rf22 =  Rf11
    Rd11 =  cosh(q.l*sqrt(abs(q.k)/(1+dE)))
    Rd12 =  sinh(q.l*sqrt(abs(q.k)/(1+dE)))/(sqrt(abs(q.k)/(1+dE)))
    Rd21 =  (sqrt(abs(q.k)/(1+dE)))*sinh(q.l*sqrt(abs(q.k)/(1+dE)))
    Rd22 =  Rd11
    if (q.k <=0 )
        println("      quad defoc")
        R=[Basic[Rd11 Rd12]; Basic[Rd21 Rd22]]
    else
        println("      quad foc")
        R=[Basic[Rf11 Rf12]; Basic[Rf21 Rf22]]
    end
    return R
end
#=
# thick quad
function mq(q)
    # horizontally  focusing quadrupole transport matrix
    # four rows
    # quad chromaticity k->k/(1+dE)
    Rf11 =   cos(q.l*sqrt(abs(q.k)/(1+dE)))
    Rf12 =  sin(q.l*sqrt(abs(q.k)/(1+dE))) / (sqrt(abs(q.k)/(1+dE)))
    Rf21 = -(sqrt(abs(q.k)/(1+dE)))*sin(q.l*sqrt(abs(q.k)/(1+dE)))
    Rf22 =  Rf11
    Rd11 =  cosh(q.l*sqrt(abs(q.k)/(1+dE)))
    Rd12 =  sinh(q.l*sqrt(abs(q.k)/(1+dE)))/(sqrt(abs(q.k)/(1+dE)))
    Rd21 =  (sqrt(abs(q.k)/(1+dE)))*sinh(q.l*sqrt(abs(q.k)/(1+dE)))
    Rd22 =  Rd11
    if (q.k <=0 )
        println("      quad defoc")
        R=[Basic[Rd11 Rd12 0 0]; Basic[Rd21 Rd22 0 0]; Basic[0 0 Rf11 Rf12]; Basic[0 0 Rf21 Rf22]]
    else
        println("      quad foc")
        R=[Basic[Rf11 Rf12 0 0]; Basic[Rf21 Rf22 0 0]; Basic[0 0 Rd11 Rd12]; Basic[0 0 Rd21 Rd22]]
    end
    return R
end
=#
function mqq(q)
    # horizontally  focusing quadrupole transport matrix
    # four rows
    # quad chromaticity k->k(1+dE)
    R=[Basic[1 0 0 0]; Basic[-(1+dE)/(q.kl) 1 0 0]; Basic[0 0 1 0];Basic[0 0 (1+dE)/(q.kl) 1]]
    return R
end
# drift 2D
function md(d)
    # four rows
#    R=[Basic[1 d.l 0 0];Basic[0 1 0 0];Basic[0 0 1 d.l];Basic[0 0 0 1]]
    println("      drift")
    R=[Basic[1 d];Basic[0 1]]
    return R
end
#=
#drfit 4D
function md(d)
    # four rows
#    R=[Basic[1 d.l 0 0];Basic[0 1 0 0];Basic[0 0 1 d.l];Basic[0 0 0 1]]
    println("      drift")
    R=[Basic[1 d 0 0];Basic[0 1 0 0];Basic[0 0 1 d];Basic[0 0 0 1]]
    return R
end
=#
### Matrix multiplication
function mult(A,B)
    p = Any[]
    ndim = 2#2D
    for i in 1:ndim, j in 1:ndim
        pp = 0
        typeof(pp)
        for k in 1:ndim
            pp = pp + dotaylor(A[i,k]*B[k,j])### order changed because push is adding column vectors
#            println("i",i," j",j," k",k," pp",pp," A",A[i,k]," B",B[k,j])
        end
        push!(p,pp) ### maybe push is adding column vectors
#        println(" stato p",p)
    end
    return transpose(reshape(p,ndim,ndim))
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
    dp  = diff(p,  dE)
    ddp = diff(dp, dE)
    f = subs(p, dE, 0) + subs(dp, dE, 0)*dE + 0.5*subs(ddp, dE, 0)*dE^2
    ftaylor = expand(f)
#=
    println("o0 : ",subs(p, dE, 0))
    println("o1 : ",subs(dp, dE, 0))
    println("o2 : ",subs(ddp, dE, 0))
=#
    return ftaylor
end
# matrix taylor expansion
function MatrixTaylor(R)
    Rarraytaylor = Any[]
    ndim = 2
    for i in 1:ndim, j in 1:ndim
        Rij = dotaylor(expand(R[i,j]))
        push!(Rarraytaylor,Rij)
    end
    Rtaylor=transpose(reshape(Rarraytaylor,ndim,ndim))
    return Rtaylor
end
