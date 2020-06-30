# thick quad matrices
function Rqf(q)
    Rf11 =  cos(q.l*sqrt(abs(q.k)/(1+dE)))
    Rf12 =  sin(q.l*sqrt(abs(q.k)/(1+dE))) / (sqrt(abs(q.k)/(1+dE)))
    Rf21 = -(sqrt(abs(q.k)/(1+dE)))*sin(q.l*sqrt(abs(q.k)/(1+dE)))
    Rf22 =  Rf11
    RF   =  [Basic[Rf11 Rf12]; Basic[Rf21 Rf22]]
    return  RF
end
function Rqd(q)
    Rd11 =  cosh(q.l*sqrt(abs(q.k)/(1+dE)))
    Rd12 =  sinh(q.l*sqrt(abs(q.k)/(1+dE)))/(sqrt(abs(q.k)/(1+dE)))
    Rd21 =  (sqrt(abs(q.k)/(1+dE)))*sinh(q.l*sqrt(abs(q.k)/(1+dE)))
    Rd22 =  Rd11
    RD   =  [Basic[Rd11 Rd12]; Basic[Rd21 Rd22]]
    return  RD
end

# quad numeric
function mq(q, verbose = 0)
    # horizontally  focusing quadrupole transport matrix
    # four rows
    # quad chromaticity k->k/(1+dE)
    if (q.k == 0 )
        if (verbose == 1 ) println("      quad k=0") end
        R=[Basic[1 q.l 0 0];Basic[0 1 0 0];Basic[0 0 1 q.l];Basic[0 0 0 1]]
    elseif (q.k < 0)
        if (verbose == 1 ) println("      quad defoc") end
        R =  [Rqd(q) zeros(Float64,2,2); zeros(Float64,2,2) Rqf(q)]
    else
        if (verbose == 1 )  println("      quad foc") end 
        R =  [Rqf(q) zeros(Float64,2,2); zeros(Float64,2,2) Rqd(q)]
    end
    return R
end

# thick quad, symbolic focusing
function mqf(q)
    # horizontally  focusing quadrupole transport matrix
    # four rows
    # quad chromaticity k->k/(1+dE)
    println("      quad foc")
    R =  [Rqf(q) zeros(Float64,2,2); zeros(Float64,2,2) Rqd(q)]
    return R
end

# thick quad, symbolic defocusing
function mqd(q)
    # horizontally  focusing quadrupole transport matrix
    # four rows
    # quad chromaticity k->k/(1+dE)
    println("      quad defoc")
    R =  [Rqd(q) zeros(Float64,2,2); zeros(Float64,2,2) Rqf(q)]
    return R
end

# drift
function md(d, verbose = 0)
    # four rows
    if (verbose == 1 ) println("      drift") end
    R=[Basic[1 d.l 0 0];Basic[0 1 0 0];Basic[0 0 1 d.l];Basic[0 0 0 1]]
#    R=[Basic[1 d];Basic[0 1]]
    return R
end

### Matrix multiplication
function mult(A,B)
    p = Any[]
    ndim = size(A,1)
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
function dotransport(mlist,verbose = 0)
    global R = mlist[1]
    for i in 1:length(mlist)-1
        if length(mlist) > 1
            #            println("enter")
            global R = mult(R,mlist[i+1]) ### mult modifies mt
        end
    end
    # Print R matrix
    if (verbose == 1 )
        for i in 1:size(R,1), j in 1:size(R,1)
            println("      R$i$j = ",R[i,j])
        end
    end
    return R
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
    ndim = size(R,1)
    for i in 1:ndim, j in 1:ndim
        Rij = dotaylor(expand(R[i,j]))
        push!(Rarraytaylor,Rij)
    end
    Rtaylor=transpose(reshape(Rarraytaylor,ndim,ndim))
    return Rtaylor
end
