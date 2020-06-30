### mpt104

function mpt104(brhofactor = 1)
    mpt104 = Any[]
    push!(mpt104,md(drift{Float64}( 0.1)))
    push!(mpt104,mq(quad{Float64}(  7/brhofactor, 0.371)))
    return mpt104
    push!(mpt104,md(drift{Float64}( 0.101)))
    push!(mpt104,mq(quad{Float64}( -7/brhofactor, 0.426)))
    push!(mpt104,md(drift{Float64}( 0.101)))
    push!(mpt104,mq(quad{Float64}(  5.759/brhofactor, 0.394)))
    push!(mpt104,md(drift{Float64}( 0.101)))
    push!(mpt104,mq(quad{Float64}( -0.543/brhofactor, 2.297)))
    push!(mpt104,md(drift{Float64}( 0.101)))
    push!(mpt104,mq(quad{Float64}(  0.761/brhofactor, 1.778)))
    push!(mpt104,md(drift{Float64}( 0.101)))
    push!(mpt104,mq(quad{Float64}( -3.388/brhofactor, 0.410)))
    push!(mpt104,md(drift{Float64}( 0.446)))
    push!(mpt104,md(drift{Float64}( 0.446)))
    push!(mpt104,mq(quad{Float64}(  3.388/brhofactor, 0.410)))
    push!(mpt104,md(drift{Float64}( 0.101)))
    push!(mpt104,mq(quad{Float64}( -0.761/brhofactor, 1.778)))
    push!(mpt104,md(drift{Float64}( 0.101)))
    push!(mpt104,mq(quad{Float64}(  0.543/brhofactor, 2.297)))
    push!(mpt104,md(drift{Float64}( 0.101)))
    push!(mpt104,mq(quad{Float64}( -5.759/brhofactor, 0.394)))
    push!(mpt104,md(drift{Float64}( 0.101)))
    push!(mpt104,mq(quad{Float64}(  7, 0.426/brhofactor)))
    push!(mpt104,md(drift{Float64}( 0.101)))
    push!(mpt104,mq(quad{Float64}( -7, 0.371/brhofactor)))
    push!(mpt104,md(drift{Float64}( 0.1 )))
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
