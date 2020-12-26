###############################################
# check if the node demand and supply capacity values are valide when there is no damage.
# Validity criteria: 1) When there is no damage, no slack is incurred.
#                    2) The flow rate of all links (acutal flow over link capacity) is no greater than 0.85.

# Nodal demand, supply capacity, and link capacity of power network are obtained from the IEEE test case (mpc object in matlab).
# link capacity of gas network is calculated according to the diameter of gas pipelines. Node demand and capacity is set to a significantly
# value and calculated by the solution on link flow according to the minCost model.
###############################################


###############################################
# import pakcages and functions
using JuMP, Gurobi
using Statistics
using CSV, DataFrames, JLD, HDF5
const GUROBI_ENV = Gurobi.Env()

cd("C:/Users/yuj5/Documents/GitHub/ICOSSAR2021")
include("./readin.jl")


###############################################
# function
# construct the min cost model to check if a valid initial feasible flow with no slack and capped flow rate exists.
function minCost(netList, networkData, interList, flowRateUB=0.85, supRate=1)
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 1))

    # decision variables
    @variable(mp, f[n in netList, k in networkData[n].inbrList]>=0)
    @variable(mp, h[interI in 1:length(interList)]>=0)
    @variable(mp, 0 <= sup[n in netList, i in networkData[n].IDList] <= networkData[n].sc[i]*supRate)

    # constraints
    # flow balance: outflow - inflow + transformed flow from other network = supply - demand
    @constraint(mp, flowbalance[n in netList, i in networkData[n].IDList],
        sum(f[n,k] for k in networkData[n].inbrList if i == k[1]) - sum(f[n,k] for k in networkData[n].inbrList if i == k[2]) +
            - sum(h[interI]*interList[interI].convRate for interI in 1:length(interList) if (interList[interI].endNet == n)&(interList[interI].endNode == i))
            == sup[n,i] - networkData[n].b[i])

    # flow capacity
    @constraint(mp, flowCap1[n in netList, k in networkData[n].inbrList],
        f[n,k] <= networkData[n].u[k]*flowRateUB)    # -networkData[n].u[k]*flowRateUB <= is removed.

    # interconnected nodes
    @constraint(mp, interConnect[interI in 1:length(interList)],
        h[interI] <= interList[interI].u)
    interDict = Dict()
    for interI in 1:length(interList)
        startN = interList[interI].startNode
        endN = interList[interI].endNode
        startNet = interList[interI].startNet
        endNet = interList[interI].endNet
        if !((endNet,endN) in keys(interDict))
            interDict[(endNet,endN)] = [interI]
        else
            push!(interDict[(endNet,endN)], interI)
        end
    end
    @constraint(mp, interCapacity[n in netList, i in networkData[n].IDList; (n,i) in keys(interDict)],
        sup[n,i] == sum(h[interI]*interList[interI].convRate for interI in 1:length(interList) if (interList[interI].endNet == n)&(interList[interI].endNode == i)))

    # objective function
    @expression(mp, supplyCost, sum(sum(networkData[n].csc[i]*sup[n,i] for i in networkData[n].IDList) for n in netList))
    @objective(mp, Min, supplyCost)

    return mp
end


function solve_and_print(netList,networkData,interList, flowRateUB=0.85, supRate=1)
    # build and solve the model
    modelMinCost = minCost(netList,networkData,interList, flowRateUB, supRate)
    optimize!(modelMinCost)

    # get results
    obj = objective_value(modelMinCost)
    f = value.(modelMinCost[:f])
    h = value.(modelMinCost[:h])
    sup = value.(modelMinCost[:sup])

    # print results
    println("\n-------------------")
    println("\nSupply at nodes\n", )
    for n in netList
        for i in networkData[n].IDList
            println(sup[n,i])
            # println(n," ",i," ",sup[n, i])
        end
    end

    println("\nFlow on arcs of each network\n")
    for n in netList
        for k in networkData[n].inbrList
            flowRate = f[n, k]/networkData[n].u[k]
            # println(n," ",k," ", f[n,k], " ", networkData[n].u[k], " ", flowRate)
            println(f[n, k])
        end
    end

    println("\nFlow on interdependent arcs\n")
    for interI in 1:length(interList)
        flowRate = h[interI]/interList[interI].u
        # println(interList[interI], " ", h[interI], " ", interList[interI].u, " ", flowRate)
        println(h[interI])
    end
end

function run_main()
    # load data from case_49 folder
    netList, networkData, interList = loadData("./data/case_49/")
    solve_and_print(netList, networkData, interList, flowRateUB=0.80, supRate==1)
end

run_main()
