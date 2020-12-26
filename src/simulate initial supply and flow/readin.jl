# read in the network data from the csv files
using CSV,DataFrames,JLD,HDF5;

# define the struct types
struct netData
    # in network
    IDList :: Array{Any,1}          # [1, 2, 3, ...]
    inbrList :: Array{Any,1}        # [(1,2), (1,3), ...]
    inbr1 :: Dict{Any,Any}          # [1]: [2, 3, ...],[2]:
    inbr2 :: Dict{Any,Any}          # [1]: [2, 3, ...],[2]:

    # demand, supply capacity, flow capacity
    b :: Dict{Any,Any}
    sc :: Dict{Any,Any}
    u :: Dict{Any,Any}
    # csc :: Dict{Any,Any}
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

    bList = Dict();
    # arc information
    brList = Dict();

    # cfaList = Dict();
    uList = Dict();
    # rtaList = Dict();
    br1List = Dict();
    br2List = Dict();
    scList = Dict();
    # cscList = Dict();
    for i in netList
        nodeList[i] = [];
        bList[i] = Dict();
        scList[i] = Dict();
        # cscList[i] = Dict();

        brList[i] = [];
        # cfaList[i] = Dict();
        uList[i] = Dict();
        br1List[i] = Dict();
        br2List[i] = Dict();
    end

    # read in the node information
    for i in 1:nnodes
        ID = nodes_data.node_id[i];
        netBelong = nodes_data.net_id[i];
        push!(nodeList[netBelong],ID);
        nodeNet[ID] = netBelong;
        bList[netBelong][ID] = nodes_data.b[i];
        scList[netBelong][ID] = nodes_data.sc[i];
        # cscList[netBelong][ID] = nodes_data.csc[i];
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
        networkData[i] = netData(nodeList[i],brList[i],br1List[i],br2List[i],bList[i],scList[i],uList[i]); #,cscList[i]);
    end

    return netList,networkData,interList;
end

function arcTrans(arcStr)
    arcStr = strip(arcStr,'(');
    arcStr = strip(arcStr,')');
    fromNodeStr,toNodeStr = split(arcStr,',');
    return (fromNodeStr,toNodeStr);
end

# load the data file given the folder address
function loadData(dataAdd)
    nodeAdd = joinpath(dataAdd,"nodes_data.csv");
    arcAdd = joinpath(dataAdd,"arcs_data.csv");
    netList,networkData,interList = readNetwork(nodeAdd,arcAdd);
    return netList,networkData,interList;
end
