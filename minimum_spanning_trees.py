
# ########## KRUSKAL BELOw
# # : Kruskal
# kruskal_g=[]
# for i in reachable:
#     for j in reachable:
#         if weight_array[i][j]!=-1:  #  I changed to be not equal to -1 : the 'not there' value
#             kruskal_g.append([i,j,weight_array[i][j]])

# # set kruskal parameters
# #result=[]
# i,e=0,0
# sorted_g = sorted(kruskal_g, key=lambda item: item[2])
# parent=[]
# rank=[]
# for node in reachable:
#     parent.append(node)
#     rank.append(0)

# # kruskal helper functions
# def kruskal_find(parent, i):
#     if parent[reachable.index(i)]==i:
#         return i
#     return kruskal_find(parent, parent[reachable.index(i)])

# def kruskal_apply_union(parent, rank, x, y):
#     xroot = kruskal_find(parent, x)
#     yroot = kruskal_find(parent, y)
#     if rank[reachable.index(xroot)] < rank[reachable.index(yroot)]:
#         parent[reachable.index(xroot)] = yroot
#     elif rank[reachable.index(xroot)] > rank[reachable.index(yroot)]:
#         parent[reachable.index(yroot)] = xroot
#     else:
#         parent[reachable.index(yroot)] = xroot
#         rank[reachable.index(xroot)] += 1

# # kruskal loop
# while e < (len(reachable) - 1):
#     u, v, w =sorted_g[i]
#     i += 1
#     x=kruskal_find(parent,u)
#     y=kruskal_find(parent,v)
#     if x != y:
#         e+=1
#         #result.append([u,v,w])
#         bush[u][v]=1
#         kruskal_apply_union(parent, rank, x, y)

# ########## KRUSKAL ABOVE

# for u, v, weight in result:
#     print("%d - %d: %d" % (u, v, weight))

################### PRIMS BELOW
# # Use Prim's algorithm to create minimum spanning tree from origin - append the bush
# prims_selected = [False]*(len(reachable))
# num_edges = 1
# # origin should be first number, so set equal to true
# prims_selected[0]=True
# # Begin Prim's algorithm
# while num_edges < (num_nodes - 1):
#     minimum=999999999999
#     x=0
#     y=0

#     for i in reachable:
#         if prims_selected[reachable.index(i)]:
#             for j in reachable:
#                 if ((weight_array[i][j] > 0) and (not prims_selected[reachable.index(j)])):
#                     if minimum > weight_array[i][j]:
#                         minimum = weight_array[i][j]
#                         x=i
#                         y=j
#     #print(str(x) + "-" + str(y) + ":" + str(weight_array[x][y]))
#     bush[x][y] = 1
#     if y == 0:
#         pass
#     else:
#         prims_selected[reachable.index(y)] = True
#     num_edges+=1
# ######################## PRIMS ABOVE







def label_function(topo_order, bush, weight_array_itter, L_link, U_link, L_node, U_node):
                 '''function to calculate all of the L and U labels for algorithm B'''
                  # determine L and U labels in forward topological ordering
                  id = 0
                   while id <= len(topo_order)-1:

                        # for first in topological order
                        # i == topo_order value & id == index
                        i = topo_order[id]
                        if id == 0:
                            # TODO: I think this can be optimized to only search for values in topo_order in position higher than what we are currently are at - the whole benefit of using topo order also check on this feature below in this same loop - CAN MAKE IT BETTER
                            for j in topo_order:
                                if bush[i][j] == 1:
                                    L_link[i][j] = weight_array_itter[i][j]
                                    U_link[i][j] = weight_array_itter[i][j]
                                    #print(i, j, U_link[i][j])

                        # for last in topological order
                        elif id == len(topo_order)-1:
                            l_val = L_node[id]
                            u_val = U_node[id]
                            for h in topo_order:
                                if bush[h][i] == 1:
                                    l_val = min(l_val, L_link[h][i])
                                    u_val = max(u_val, U_link[h][i])
                            L_node[id] = l_val
                            U_node[id] = u_val

                        # for all others
                        else:

                            l_val = L_node[id]
                            u_val = U_node[id]
                            # is this the correct slicing?
                            for h in topo_order[:id+1]:

                                if bush[h][i] == 1:
                                    l_val = min(l_val, L_link[h][i])
                                    u_val = max(u_val, U_link[h][i])

                            # OLD METHOD BELOW
                            # L_node[id] = l_val
                            # U_node[id] = u_val
                            # for j in topo_order[id:]:  # is this the correct slicing
                            #     if bush[i][j] == 1:
                            #         L_link[i][j] = min(L_link[i][j], L_node[id]+weight_array_itter[i][j])
                            #         U_link[i][j] = max(U_link[i][j], U_node[id]+weight_array_itter[i][j])
                            # OLD METHOD ABOVE

                            # begin BUG fix (replace lines above with below)
                            #  L and U Labels being calculated incorrectly
                            # potentially because multiple nodes may have zero inputs (compared to test data where only 1 node had zero nodes going into itself)
                            if l_val == np.inf:
                                for j in topo_order[id:]:  # is this correct slicing
                                    if bush[i][j] == 1:
                                        L_link[i][j] = weight_array_itter[i][j]
                                        U_link[i][j] = weight_array_itter[i][j]
                                # : are below two lines necessary? Are they just wrong?
                                L_node[id] = 0
                                U_node[id] = 0

                            else:
                                L_node[id] = l_val
                                U_node[id] = u_val
                                # is this the correct slicing
                                for j in topo_order[id:]:
                                    if bush[i][j] == 1:
                                        L_link[i][j] = min(L_link[i][j],
                                                           L_node[id]+weight_array_itter[i][j])
                                        U_link[i][j] = max(U_link[i][j],
                                                           U_node[id]+weight_array_itter[i][j])

                            ## END BUG FIX

                        id += 1

                    # remove infs for clarity
                    U_link[np.isinf(U_link)] = 0
                    L_link[np.isinf(L_link)] = 0

                    return L_node, U_node, L_link, U_link
