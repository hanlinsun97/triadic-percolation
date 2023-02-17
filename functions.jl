
 # If you use this code, please cite:
 # Hanlin Sun, Filippo Radicchi, Jürgen Kurths and Ginestra Bianconi
 # "The dynamic nature of percolation on networks with triadic interactions"
 # arXiv:2204.13067


using SpecialFunctions
using Random
using LinearAlgebra

# All the functions used in the simulation are defined here.

function sf_network_hard_constraint(N, gamma, m, cp, cn)

    # Generate a uncorrelated degree sequence with N nodes, given the power-law degree distribution
    # with minimum degree m and exponent gamma

    # This is a configuration model, the constraints on degrees are hard.

    max_degree = sqrt(N)    #  Structural cutoff

    prob = collect(m:max_degree) .^ (-gamma)
    prob = prob ./ sum(prob)

    degree_seq = zeros(Int64, N)

    for i = 1:N
        comp = [0; cumsum(prob)] .< rand()
        j = 1
        while (comp[j] == 1)
            j += 1
        end
        degree_seq[i] = j - 1 + m - 1
    end

    while iseven(sum(degree_seq)) == 0
        # The sum of degree of all nodes of a network must be even.
        # If not, resample
        for i = 1:N
            comp = [0; cumsum(prob)] .< rand()
            j = 1
            while (comp[j] == 1)
                j += 1
            end
            degree_seq[i] = j - 1 + m - 1
        end
    end

    # Generate a network using the input degree sequence

    # Generate two dictionaries to store the network
    # dict_nodes_stru: nodeID -> [edgeID of the neighbor of this node]
    # dict_edges_stru: edgeID -> [nodeID of the neighbor of this edge]

    dict_nodes_stru = Dict{Int64,Vector{Int64}}()
    dict_edges_stru = Dict{Int64,Vector{Int64}}()

    M = sum(degree_seq) / 2    # number of links
    N = length(degree_seq)     # number of nodes

    for i = 1:N
        dict_nodes_stru[i] = Vector{Int64}()
    end

    for i = N+1:N+M
        dict_edges_stru[i] = Vector{Int64}()
    end

    edge_counter = degree_seq
    alive_node = collect(1:N)   # nodes that still have free studs
    edgeID = N + 1

    while sum(edge_counter) > 0

        i = rand(1:length(alive_node))
        node1 = alive_node[i]
        edge_counter[node1] -= 1    # number of free studs remains
        filter!(x -> edge_counter[x] !== 0, alive_node)

        j = rand(1:length(alive_node))
        node2 = alive_node[j]
        edge_counter[node2] -= 1
        filter!(x -> edge_counter[x] !== 0, alive_node)

        push!(dict_nodes_stru[node1], edgeID)
        push!(dict_nodes_stru[node2], edgeID)
        push!(dict_edges_stru[edgeID], node1)
        push!(dict_edges_stru[edgeID], node2)

        edgeID += 1
    end

    # Generate the Poissonian regulatory network

    dict_edges_pos = Dict{Int64,Vector{Int64}}()
    dict_edges_neg = Dict{Int64,Vector{Int64}}()

    N = length(dict_nodes_stru)
    M = length(dict_edges_stru)

    for i = N+1:N+M
        dict_edges_pos[i] = Int64[]
        dict_edges_neg[i] = Int64[]
    end

    for edge = N+1:N+M
        for node = 1:N
            if !(node in dict_edges_stru[edge])
                if rand() < cp / N
                    push!(dict_edges_pos[edge], node)
                end
            end
        end
    end

    for edge = N+1:N+M
        for node = 1:N
            if !(node in dict_edges_stru[edge]) & !(node in dict_edges_pos[edge])
                if rand() < cn / N
                    push!(dict_edges_neg[edge], node)
                end
            end
        end
    end
    return dict_nodes_stru, dict_edges_stru, dict_edges_pos, dict_edges_neg
end

function sf_network_soft_constraint(N, gamma, m, cp, cn)

    # Structural Scale-free network with size N, exponent gamma and minimum degree m.
    # Positive and negative regulatory networks are Poisson networks with average degree cp and cn.
    # The constraints on node degree are soft.
    dict_nodes_stru = Dict{Int64,Vector{Int64}}()
    dict_edges_stru = Dict{Int64,Vector{Int64}}()
    dict_edges_pos = Dict{Int64,Vector{Int64}}()
    dict_edges_neg = Dict{Int64,Vector{Int64}}()

    upper_bound = N * (zeta(gamma - 1, m) / zeta(gamma, m)) * 2  # preallocation, max number of edges

    for i = 1:N
        dict_nodes_stru[i] = Int64[]
    end

    for j = N+1:N+upper_bound
        dict_edges_stru[j] = Int64[]
    end
    norm = 0
    cutoff = sqrt(N)
    kx1 = zeros(Float64, N)
    for i = 1:N
        kx1[i] = (rand() * m^(1 - gamma))^(-1 / (gamma - 1))
        while (kx1[i] > cutoff)
            kx1[i] = (rand() * m^(1 - gamma))^(-1 / (gamma - 1))
        end
        norm += kx1[i]
    end

    edgeID = N + 1
    for i = 1:N
        for j = i+1:N
            if rand() < (kx1[i] * kx1[j]) / norm
                push!(dict_nodes_stru[i], edgeID)
                push!(dict_nodes_stru[j], edgeID)
                push!(dict_edges_stru[edgeID], i)
                push!(dict_edges_stru[edgeID], j)
                edgeID += 1
            end
        end
    end

    max_edgeID = edgeID - 1
    for i = max_edgeID+1:upper_bound+N
        delete!(dict_edges_stru, i)
    end
    M = max_edgeID - N  # Num of edges

    for i = N+1:N+M
        dict_edges_pos[i] = Int64[]
        dict_edges_neg[i] = Int64[]
    end

    for edge = N+1:N+M
        for node = 1:N
            if !(node in dict_edges_stru[edge])
                if rand() < cp / N
                    push!(dict_edges_pos[edge], node)
                end
            end
        end
    end

    for edge = N+1:N+M
        for node = 1:N
            if !(node in dict_edges_stru[edge]) & !(node in dict_edges_pos[edge])
                if rand() < cn / N
                    push!(dict_edges_neg[edge], node)
                end
            end
        end
    end
    return dict_nodes_stru, dict_edges_stru, dict_edges_pos, dict_edges_neg
end



function poisson_network(N, c, cp, cn)

    # Structural Poisson network with size N and average degree c.
    # Positive and negative regulatory networks are Poisson networks with average degree cp and cn.

    dict_nodes_stru = Dict{Int64,Vector{Int64}}()
    dict_edges_stru = Dict{Int64,Vector{Int64}}()
    dict_edges_pos = Dict{Int64,Vector{Int64}}()
    dict_edges_neg = Dict{Int64,Vector{Int64}}()

    upper_bound = N * c * 2  # preallocation, max number of edges

    for i = 1:N
        dict_nodes_stru[i] = Int64[]
    end

    for j = N+1:N+upper_bound
        dict_edges_stru[j] = Int64[]
    end

    edgeID = N + 1
    for i = 1:N
        for j = i+1:N
            if rand() < c / N
                push!(dict_nodes_stru[i], edgeID)
                push!(dict_nodes_stru[j], edgeID)
                push!(dict_edges_stru[edgeID], i)
                push!(dict_edges_stru[edgeID], j)
                edgeID += 1
            end
        end
    end

    max_edgeID = edgeID - 1
    for i = max_edgeID+1:upper_bound+N
        delete!(dict_edges_stru, i)
    end

    M = max_edgeID - N  # Num of edges

    for i = N+1:N+M
        dict_edges_pos[i] = Int64[]
        dict_edges_neg[i] = Int64[]
    end

    for edge = N+1:N+M
        for node = 1:N
            if !(node in dict_edges_stru[edge])
                if rand() < cp / N
                    push!(dict_edges_pos[edge], node)
                end
            end
        end
    end

    for edge = N+1:N+M
        for node = 1:N
            if !(node in dict_edges_stru[edge]) & !(node in dict_edges_pos[edge])
                if rand() < cn / N
                    push!(dict_edges_neg[edge], node)
                end
            end
        end
    end
    return dict_nodes_stru, dict_edges_stru, dict_edges_pos, dict_edges_neg
end


function bfs(node, dict_nodes, dict_edges, counted, retained, clusterID, ID_list)

    # Breadth first search algorithm
    # Given the network structure, return the size of component that `node` belongs to
    # and assign this component an ID called `clusterID`

    # counted: 0 not counted, 1 counted
    # retained: 0 not retained, 1 retained

    queue_edges = Int64[]
    queue_nodes = Int64[]

    push!(queue_nodes, node)

    counted[node] = 1
    ID_list[node] = clusterID
    cluster_size = 1

    while !isempty(queue_nodes)
        node = popfirst!(queue_nodes)

        for edge in dict_nodes[node]
            if (counted[edge] == 0) & (retained[edge] == 1)
                push!(queue_edges, edge)
                counted[edge] = 1
                ID_list[edge] = clusterID
            end
        end

        while !isempty(queue_edges)
            edge = popfirst!(queue_edges)
            for node in dict_edges[edge]
                if (counted[node] == 0) & (retained[node] == 1)
                    cluster_size += 1
                    push!(queue_nodes, node)
                    counted[node] = 1
                    ID_list[node] = clusterID
                end
            end
        end
    end
    return cluster_size, counted, ID_list, clusterID
end


function gc(dict_nodes, dict_edges, retained)

    # Using the BFS algorithm, given the network structure, return the size of the largest component and its ID

    num_nodes = length(dict_nodes)
    num_edges = length(dict_edges)
    counted = zeros(Int64, num_nodes + num_edges)
    ID_list = zeros(Int64, num_nodes + num_edges)
    max_cluster_size = 0
    max_clusterID = 0
    clusterID = 0
    for i = 1:num_nodes
        if (counted[i] == 0) & (retained[i] == 1)
            clusterID += 1
            cluster_size, counted, ID_list, clusterID = bfs(i, dict_nodes, dict_edges, counted, retained, clusterID, ID_list)
            if cluster_size > max_cluster_size
                max_cluster_size = cluster_size
                max_clusterID = clusterID
            end
        else
            continue
        end
    end
    return max_cluster_size, max_clusterID, ID_list
end


function sf_theory(degree_seq, cp, cn)

    # The theoretical orbit diagram of p versus R, using the generating function method.

    k = degree_seq
    N = length(k)
    total_time = 300    # Maximum time step
    n = 150     # Take the last n samples from the time series

    R_collection = Int64[]
    p_collection = Int64[]
    
    for p = 1:-0.01:0
        SN2 = 0.99
        pL = 0.1
        R = zeros(total_time)
        timeR = zeros(total_time)

        for i = 1:total_time

            for nrun = 1:30
                SN2 = (1 .- sum(k .* (1 - pL * SN2) .^ (k .- 1)) / sum(k))
            end

            SN1 = (1 .- sum((1 - pL * SN2) .^ k) / N)
            pL = p * abs(exp(-cn * SN1) * (1 - exp(-cp * SN1)))
            R[i] = SN1
            timeR[i] = i
            # @show SN1, i
        end

        R_collection = [R_collection; R[end-n+1:end]]   # Take the last n samples from the time series
        p_collection = [p_collection; p * ones(n)]
    end
    return R_collection, p_collection
end


function poisson_theory(c, cp, cn)

    R_collection = []
    p_collection = []

    total_time = 100       # Maximum time step
    n = 50      # Take the last n samples from the time series
    
    R = zeros(total_time)
    timeR = zeros(total_time)
    SL2 = 1
    SN2 = 1
    pL = 0.9

    for p = 1:-0.01:0

        for i = 1:total_time
            for nrun = 1:100
                SN2 = (1 - exp(-c * SL2))
                SL2 = pL * SN2
            end
            pL = p * abs(exp(-cn * SN2) * (1 - exp(-cp * SN2)))
            R[i] = SN2
            timeR[i] = i
        end
        R_collection = [R_collection; R[end-n+1:end]]   # Take the last n samples from the time series
        p_collection = [p_collection; p * ones(n)]
    end
    return p_collection, R_collection
end


function simulation(dict_nodes_stru, dict_edges_stru, dict_edges_pos, dict_edges_neg, p, retained, Tmax)

    #  Numerical simulation on the given network
    #  Retained: initial condition of edges. 1: retained, 0: damaged
    #  Tmax: duration of the time series

    N = length(dict_nodes_stru)
    
    active_node = ones(N)

    # Finding the giant component and deactivate the nodes and links which are not in the giant component.
    
    # Collect the time series
    cluster_size_list = zeros(Tmax)
    T_list = zeros(Tmax)

    for t = 1:Tmax

        # Find the giant component
        max_cluster_size, max_clusterID, ID_list = gc(dict_nodes_stru, dict_edges_stru, retained)

        for i = 1:N
            if ID_list[i] != max_clusterID
                active_node[i] = 0      # Only nodes in the giant component are active.
            else
                active_node[i] = 1
            end
        end
        cluster_size_list[t] = (max_cluster_size / N)
        T_list[t] = t
      
        # Update the pL
        for edge in keys(dict_edges_stru)

            pos_regulation = dict_edges_pos[edge]
            neg_regulation = dict_edges_neg[edge]
            
            sum_pos = 0
            sum_neg = 0
          

            for pos_node in pos_regulation
                sum_pos += active_node[pos_node]
            end

            for neg_node in neg_regulation
                sum_neg += active_node[neg_node]
            end

            if (sum_pos > 0) & (sum_neg == 0) & (rand() < p) 
                retained[edge] = 1
            else
                retained[edge] = 0
            end
        end
    end
    return cluster_size_list, T_list, retained
end


function amp_theory_p(c, cp, cn)
    
    # Calculation of the amplitude of the period-2 oscillation, as a function of p.

    amp_collection = []
    p_collection = []

    total_time = 2000
    R = zeros(total_time)
    timeR = zeros(total_time)
    for p = 0.5:0.001:0.8
        @show p
        SL2=1
        SN2=1
        pL=p
        for i=1:total_time
            for nrun=1:1000
                SN2=(1-exp(-c*SL2))
                SL2=pL*SN2
            end
            pL= p * abs(exp(-cn*SN2)*(1-exp(-cp*SN2)))
            R[i]=SN2;
            timeR[i]=i;
        end
        push!(amp_collection, abs(R[end]-R[end-1]))
        push!(p_collection, p)
    end
    return p_collection, amp_collection
end

function poisson_only_negative(c, cn)

    # Theoretical orbit diagram on Poisson network when only negative regulatory interactions are present
    # (In this situation, an edge is active if it does not connect to any active negative regulators)

    R_collection = []
    p_collection = []

    N = 500
    total_time = 2000
    R = zeros(total_time)
    timeR = zeros(total_time)
    for p = 0:0.01:1
        @show p
        SL2 = 1
        SN2 = 1
        pL = p
        for i = 1:total_time
            for nrun = 1:1000
                SN2 = (1 - exp(-c * SL2))
                SL2 = pL * SN2
            end
            pL = p * abs(exp(-cn * SN2))
            R[i] = SN2
            timeR[i] = i
        end
        R_collection = [R_collection; R[end-N+1:end]]
        p_collection = [p_collection; p * ones(N)]
    end
    return p_collection, R_collection
end
