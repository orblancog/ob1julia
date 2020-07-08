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
# twiss formulation
function twisspropagation(R,bx0,by0,ax0,ay0,gx0,gy0)
    Rtaylor = R
    println("    Taylor expansion to second order of the matrix R")
    ### taylor expansion already done during multiplication
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
    exprbetx0=expand(subs(subs(subs(eqbetaxdE0, alfax0, ax0), gamax0, gx0),betax0, bx0))
    exprbetx1=expand(subs(subs(subs(eqbetaxdE1, alfax0, ax0), gamax0, gx0),betax0, bx0))
    exprbetx2=expand(subs(subs(subs(eqbetaxdE2, alfax0, ax0), gamax0, gx0),betax0, bx0))
    expralfx0=expand(subs(subs(subs(eqalfaxdE0, alfax0, ax0), gamax0, gx0),betax0, bx0))
    expralfx1=expand(subs(subs(subs(eqalfaxdE1, alfax0, ax0), gamax0, gx0),betax0, bx0))
    expralfx2=expand(subs(subs(subs(eqalfaxdE2, alfax0, ax0), gamax0, gx0),betax0, bx0))
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
exprbety0=expand(subs(subs(subs(eqbetaydE0, alfay0, ay0), gamay0, gy0),betay0, by0))
exprbety1=expand(subs(subs(subs(eqbetaydE1, alfay0, ay0), gamay0, gy0),betay0, by0))
exprbety2=expand(subs(subs(subs(eqbetaydE2, alfay0, ay0), gamay0, gy0),betay0, by0))
expralfy0=expand(subs(subs(subs(eqalfaydE0, alfay0, ay0), gamay0, gy0),betay0, by0))
expralfy1=expand(subs(subs(subs(eqalfaydE1, alfay0, ay0), gamay0, gy0),betay0, by0))
expralfy2=expand(subs(subs(subs(eqalfaydE2, alfay0, ay0), gamay0, gy0),betay0, by0))
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
return exprbety0,exprbety1,exprbety2

end
