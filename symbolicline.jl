## symbolic line


Nq = 10 # number of quadrupoles
Nl = Nq # number of drifts
#println("    Using ",Nq," quads and "Nl," drifts")

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
listq = Any[]
listd = Any[]
for i in 1:Nq
    push!(listq,quad(kk[i],ll[i])) #quads
    push!(listd,drift(dd[i]))
end
println("    elements created.")


### create matrices
println("    Creating matrix representation of the transport line...")
mlist = Any[]


mlist4 = Any[]
#push!(mlist4,md(dstar))
#push!(mlist4,mqf(listq[1]))
# push!(mlist4,md(listd[1]))
# push!(mlist4,mqd(listq[2]))
# push!(mlist4,md(listd[2]))
# push!(mlist4,mqf(listq[3]))
# push!(mlist4,md(listd[3]))
# push!(mlist4,mqd(listq[4]))
# push!(mlist4,md(dstar))



# #mdstar = md(dstar)
# mdstar = md(lstar)
# push!(mlist,mdstar) # last element first
# #push!(mlist, mq(quad(kk[Nq],ll[Nq])))
# push!(mlist, mqq(qquad(kl[Nq])))
# for i in Nl:-1:1
# #    push!(mlist, md(drift(dd[i])))
#     push!(mlist, md(1.0))
# #    push!(mlist, mq(quad(kk[i],ll[i])))
#     push!(mlist, mqq(qquad(kl[i])))
# end
# push!(mlist,mdstar) # first elemen last 
# println("    ",length(mlist)," matrices created.")
# println("    Transport line created. Adios !")

#mlist2 = Any[]
#push!(mlist2,mq(quad(kk[1],ll[1])))
#push!(mlist2,md(1))
#push!(mlist2,mq(quad(-7,1)))
#push!(mlist2,md(1))
#push!(mlist2,mq(quad(7,1)))
#push!(mlist2,mdstar)
# #push!(mlist2,mq(quad(kk[2],ll[2])))
# #push!(mlist2,mdstar)
# println("length",length(mlist2))
# println(mlist2)
# exit()
