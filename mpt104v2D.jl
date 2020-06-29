### mpt104

function mpt104v2D()
    mpt104 = Any[]
    push!(mpt104,md(0.1))
    push!(mpt104,mq(quad( 7,0.371)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad(-7,0.426)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad( 5.759,0.394)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad(-0.543,2.297)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad( 0.761,1.778)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad(-3.388,0.410)))
    push!(mpt104,md(0.446))
    push!(mpt104,md(0.446))
    push!(mpt104,mq(quad( 3.388,0.410)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad(-0.761,1.778)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad( 0.543,2.297)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad(-5.759,0.394)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad( 7,0.426)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad(-7,0.185405*2)))
    push!(mpt104,md(0.1))
    return mpt104
end

#=
## careful, test only
function mpt104()
    mpt104 = Any[]
    push!(mpt104,md(0.1))
    push!(mpt104,mq(quad(-sqrt(7),0.371)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad( sqrt(7),0.426)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad(-sqrt(5.759),0.394)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad( sqrt(0.543),2.297)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad(-sqrt(0.761),1.778)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad( sqrt(3.388),0.410)))
    push!(mpt104,md(0.446))
    push!(mpt104,md(0.446))
    push!(mpt104,mq(quad(-sqrt(3.388),0.410)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad( sqrt(0.761),1.778)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad(-0.543,2.297)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad( sqrt(5.759),0.394)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad(-sqrt(7),0.426)))
    push!(mpt104,md(0.101))
    push!(mpt104,mq(quad( sqrt(7),0.371)))
    push!(mpt104,md(0.1))
    return mpt104
end
=#