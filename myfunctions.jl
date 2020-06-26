
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
# thick quad
function mq(q)
    # horizontally  focusing quadrupole transport matrix
    # four rows
    # quad chromaticity k->k/(1+dE)
    Rf11=cos(q.l*q.k/(1+dE))
    Rf12=sin(q.l*q.k/(1+dE))/(q.k/(1+dE))
    Rf21=-(q.k/(1+dE))*sin(q.l*q.k/(1+dE))
    Rf22=cos(q.l*q.k/(1+dE))
    Rd11=cosh(q.l*q.k/(1+dE))
    Rd12=sinh(q.l*q.k/(1+dE))/(q.k/(1+dE))
    Rd21=(q.k/(1+dE))*sinh(q.l*q.k/(1+dE))
    Rd22=cosh(q.l*q.k/(1+dE))
#!    if (q.k <=0 )
#        println("defoc")
#        R=[Basic[Rd11 Rd12 0 0]; Basic[Rd21 Rd22 0 0]; Basic[0 0 Rf11 Rf12]; Basic[0 0 Rf21 Rf22]]
#    else
#        println("foc")
        R=[Basic[Rf11 Rf12 0 0]; Basic[Rf21 Rf22 0 0]; Basic[0 0 Rd11 Rd12]; Basic[0 0 Rd21 Rd22]]
#    end
    return R
end
function mqq(q)
    # horizontally  focusing quadrupole transport matrix
    # four rows
    # quad chromaticity k->k(1+dE)
    R=[Basic[1 0 0 0]; Basic[-(1+dE)/(q.kl) 1 0 0]; Basic[0 0 1 0];Basic[0 0 (1+dE)/(q.kl) 1]]
    return R
end
function md(d)
    # four rows
#    R=[Basic[1 d.l 0 0];Basic[0 1 0 0];Basic[0 0 1 d.l];Basic[0 0 0 1]]
    R=[Basic[1 d 0 0];Basic[0 1 0 0];Basic[0 0 1 d];Basic[0 0 0 1]]
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
### truncate expansion
function dotruncate(p)
    pexpnd = expand(p)
    ptrunc = 0
    for i in 0:1:2 ## zeroth to second order expansion
        ptrunc = ptrunc + coeff(pexpnd,dE,Basic(i))*dE^i
    end
    return expand(ptrunc)
end

