# define the struct types
struct netData
    # in network
    IDList :: Array{Any,1}          # [1, 2, 3, ...]
    inbrList :: Array{Any,1}        # [(1,2), (1,3), ...]
    inbr1 :: Dict{Any,Any}          # [1]: [2, 3, ...],[2]:
    inbr2 :: Dict{Any,Any}          # [1]: [2, 3, ...],[2]:

    # restoration time
    rtn :: Dict{Any,Any}
    rta :: Dict{Any,Any}

    # costs
    chn :: Dict{Any,Any}
    cha :: Dict{Any,Any}
    csn :: Dict{Any,Any}
    crn :: Dict{Any,Any}
    cra :: Dict{Any,Any}

    # demand, supply capacity, flow capacity
    b :: Dict{Any,Any}
    sc :: Dict{Any,Any}
    u :: Dict{Any,Any}
    csc :: Dict{Any,Any}
end

# data structure for the inter network arcs
struct interData
    startNet :: Int64
    endNet :: Int64
    startNode :: Any
    endNode :: Any
    convRate :: Float64
    u :: Float64
end

# data structure for the scenarios
struct scenarioData
    dNodes :: Array{Any,1}
    dArcs :: Array{Any,1}
end
