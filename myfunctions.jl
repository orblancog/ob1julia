
### Matrix declaration
function mq(q)
    # horizontally  focusing quadrupole transport matrix
    # four rows
    # quad chromaticity k->k(1+dE)
    R=[Basic[1 0 0 0]; Basic[-(1+dE)/(q.k*q.l) 1 0 0]; Basic[0 0 1 0];Basic[0 0 (1+dE)/(q.k*q.l) 1]]
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
#    print("create")
    fstart   = p#expand(p)
    f = 0
    for i in 0:1:2 ## zeroth to second order expansion
        f = f + coeff(fstart,dE,Basic(i))*dE^i
    end
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

