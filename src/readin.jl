# read in the network data from the csv files
using CSV,DataFrames,JLD,HDF5;

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



# define functions

function readNetwork(nodeAdd,arcAdd)
    nodes_data = CSV.read(nodeAdd);
    arcs_data = CSV.read(arcAdd);
    nnodes,mnodes = size(nodes_data);
    narcs,marcs = size(arcs_data);

    # initiate the data structure
    # the keys to the dictionaries are the network number
    netList = unique(nodes_data.net_id);
    nodeNet = Dict();
    arcNet = Dict();
    # node information
    nodeList = Dict();
    chnList = Dict();
    crnList = Dict();
    csnList = Dict();
    rtnList = Dict();
    bList = Dict();
    # arc information
    brList = Dict();
    chaList = Dict();
    craList = Dict();
    cfaList = Dict();
    uList = Dict();
    rtaList = Dict();
    br1List = Dict();
    br2List = Dict();
    scList = Dict();
    cscList = Dict();
    for i in netList
        nodeList[i] = [];
        chnList[i] = Dict();
        crnList[i] = Dict();
        csnList[i] = Dict();
        rtnList[i] = Dict();
        bList[i] = Dict();
        scList[i] = Dict();
        cscList[i] = Dict();

        brList[i] = [];
        chaList[i] = Dict();
        craList[i] = Dict();
        cfaList[i] = Dict();
        uList[i] = Dict();
        rtaList[i] = Dict();
        br1List[i] = Dict();
        br2List[i] = Dict();
    end

    # read in the node information
    for i in 1:nnodes
        ID = nodes_data.node_id[i];
        netBelong = nodes_data.net_id[i];
        push!(nodeList[netBelong],ID);
        nodeNet[ID] = netBelong;
        chnList[netBelong][ID] = nodes_data.chn[i];
        crnList[netBelong][ID] = nodes_data.crn[i];
        csnList[netBelong][ID] = nodes_data.csn[i];
        rtnList[netBelong][ID] = nodes_data.rn[i];
        bList[netBelong][ID] = nodes_data.b[i];
        scList[netBelong][ID] = nodes_data.sc[i];
        cscList[netBelong][ID] = nodes_data.csc[i];
    end

    # read in the arc information
    interList = [];
    for a in 1:narcs
        fromNode = arcs_data.start_node[a];
        toNode = arcs_data.end_node[a];
        arcID = (fromNode,toNode);
        convRate = arcs_data.conv_rate[a];
        if nodeNet[fromNode] == nodeNet[toNode]
            # if it is within some network
            netBelong = nodeNet[fromNode];
            push!(brList[netBelong],arcID);
            arcNet[arcID] = netBelong;
            chaList[netBelong][arcID] = arcs_data.cha[a];
            craList[netBelong][arcID] = arcs_data.cra[a];
            rtaList[netBelong][arcID] = arcs_data.ra[a];
            uList[netBelong][arcID] = arcs_data.u[a];
            if fromNode in keys(br1List[netBelong])
                push!(br1List[netBelong][fromNode],arcID);
            else
                br1List[netBelong][fromNode] = [arcID];
            end
            if toNode in keys(br2List[netBelong])
                push!(br2List[netBelong][toNode],arcID);
            else
                br2List[netBelong][toNode] = [arcID];
            end
        else
            # if it is inter network
            arcInfo = interData(nodeNet[fromNode],nodeNet[toNode],fromNode,toNode,convRate,arcs_data.u[a]);
            push!(interList,arcInfo);
        end
    end

    # create networkData: a list of networks
    networkData = Dict();
    for i in netList
        networkData[i] = netData(nodeList[i],brList[i],br1List[i],br2List[i],rtnList[i],rtaList[i],
            chnList[i],chaList[i],csnList[i],crnList[i],craList[i],bList[i],scList[i],uList[i],cscList[i]);
    end

    return netList,networkData,interList;
end

function arcTrans(arcStr)
    arcStr = strip(arcStr,'(');
    arcStr = strip(arcStr,')');
    fromNodeStr,toNodeStr = split(arcStr,',');
    return (fromNodeStr,toNodeStr);
end

# read the data from a .jld file
function loadDataJLD(caseDataAdd)
    dataRaw = load(caseDataAdd);
    netList = dataRaw["netList"];
    networkData = dataRaw["networkData"];
    interList = dataRaw["interList"];
    return netList,networkData,interList;
end

# load the data file given the folder address
function loadData(dataAdd)
    nodeAdd = joinpath(dataAdd,"nodes_data.csv");
    arcAdd = joinpath(dataAdd,"arcs_data.csv");
    dNodeAdd = joinpath(dataAdd,"is_damage_nodes.csv");
    dArcAdd = joinpath(dataAdd,"is_damage_arcs.csv");
    netList,networkData,interList = readNetwork(nodeAdd,arcAdd);
    return netList,networkData,interList;
end
