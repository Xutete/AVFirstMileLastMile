# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 06:30:50 2020

@author: Pramesh Kumar
"""
import math, time
from gurobipy import *
from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

################################################################################################+

class Zone:
    def __init__(self, _tmpIn):
        self.lat = float(_tmpIn[0])
        self.long = float(_tmpIn[1])
        self.drivers = 10
        self.dest = []
        self.area = 2408689
        self.wait = 0
        
        
class Node:
    def __init__(self, _tmpIn):
        self.lat = float(_tmpIn[0])
        self.long = float(_tmpIn[1])
        self.type = _tmpIn[2]
        self.outLinks = []
        self.inLinks = []
        self.label = 0
        self.pred = ""
        self.name = ""
            
class Link:
    def __init__(self, _tmpIn, _type):
        self.fromNode = _tmpIn[0]
        self.toNode = _tmpIn[1]
        self.dist = round(float(_tmpIn[2]),2) # in miles
        self.type = _type
        
        if _type == 'zoneAccess':
            self.time =  round(transitFare*60/VOT) # in minutes
        elif _type == 'zoneEgress':
            self.time = 0
        elif _type == 'road':
            self.time =  round(float(_tmpIn[3])/0.6+ float(_tmpIn[3])*fuelCost*60/(VOT* 0.6)) # in minutes 0.6 is from 0.1 hours as the units
        else:
            self.time = round(float(_tmpIn[3]))
        self.wait = 0
        self.fuelCost = 0
        self.passengers = []
        self.active = 1
        self.lineId = ""
        
        
class Line:
    def __init__(self, _tmpIn):
        self.lineId = _tmpIn[0]
        self.lineName = _tmpIn[1]
        self.trips = [] # Associated trips
        self.stops = []
        self.links = []
        self.freq = 60 
        self.active = 1
        self.cost = 1000

class Passenger:
    def __init__(self, _tmpIn):
        self.origin = _tmpIn[0]
        self.dest = _tmpIn[1]
        self.path = []



class Demand:
    def __init__(self, _tmpIn):
        self.fromZone = _tmpIn[0]
        self.toZone = _tmpIn[1]
        self.demand = float(_tmpIn[2])/10

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    c = 2 * math.asin(math.sqrt(a))
    mi = 3959 * c
    return mi


################################################################################################


def readNodes():
    """
    Read zones, road nodes, and transit stopID
    """
        
    # Reading road nodes
    inFile = open(loc+"ft_input_zones.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        if tmpIn[0] not in nodeSet:
            nodeSet[tmpIn[0]] = Node([tmpIn[1], tmpIn[2], "road"])
            zoneSet[tmpIn[0]] =Zone(tmpIn)
        else:
            print(tmpIn[2], " roadNode already present as ", nodeSet[tmpIn[2]].type)
    inFile.close()
        
    # Reading transit nodes
    inFile = open(loc+"ft_input_stops.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        if tmpIn[0] not in nodeSet:
            nodeSet[tmpIn[0]] = Node([tmpIn[3], tmpIn[4], 'stop'])
        else:
            print(tmpIn[0], " stop already present as ", nodeSet[tmpIn[0]].type)
    inFile.close()


def readLinks():
    """
    Read access, mode transfer, road, and
    transit transfer links
    """
    # Reading access links
    inFile = open(loc+"ft_input_accessLinks.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        linkId =(tmpIn[0], tmpIn[1])
        if linkId not in linkSet:
            linkSet[linkId] = Link(tmpIn, "zoneAccess")
            if linkId not in nodeSet[tmpIn[0]].outLinks:
                nodeSet[tmpIn[0]].outLinks.append(linkId)
            if linkId not in nodeSet[tmpIn[1]].inLinks:
                nodeSet[tmpIn[1]].inLinks.append(linkId)
        else:
            print(linkId, " stop already present as ", linkSet[linkId].type)
        linkId =(tmpIn[1], tmpIn[0])
        if linkId not in linkSet:
            linkSet[linkId] = Link([tmpIn[1], tmpIn[0], tmpIn[2], tmpIn[3]], "zoneEgress")
            if linkId not in nodeSet[tmpIn[1]].outLinks:
                nodeSet[tmpIn[1]].outLinks.append(linkId)
            if linkId not in nodeSet[tmpIn[0]].inLinks:
                nodeSet[tmpIn[0]].inLinks.append(linkId)
        else:
            print(linkId, " stop already present as ", linkSet[linkId].type)
    inFile.close()


    # Reading road links
    inFile = open(loc+"network.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        linkId = (tmpIn[0], tmpIn[1])
        if linkId not in linkSet:
            linkSet[linkId] = Link([tmpIn[0], tmpIn[1], tmpIn[3], tmpIn[4]], "road")
            if linkId not in nodeSet[tmpIn[0]].outLinks:
                nodeSet[tmpIn[0]].outLinks.append(linkId)
            if linkId not in nodeSet[tmpIn[1]].inLinks:
                nodeSet[tmpIn[1]].inLinks.append(linkId)
    inFile.close()
    
def readLines():
    # Reading transit stops
    inFile = open(loc + "ft_input_routes.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        lineSet[tmpIn[0]] = Line([tmpIn[0], tmpIn[2]])
    inFile.close()
    
def readtrips():
    # Reading transit stops
    inFile = open(loc + "ft_input_trips.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        lineSet[tmpIn[1]].trips.append(tmpIn[0])
    inFile.close()
    
    
def readTransitLinks():
    inFile = open(loc+"ft_input_stopTimes.dat")
    tmpIn = inFile.readline().strip().split("\t")
    prevNodeId = ""
    for x in inFile:
        tmpIn = x.strip().split("\t")
        tripId = tmpIn[0]
        routeId = [k for k in lineSet if tripId in lineSet[k].trips][0]
        nodeId = (tmpIn[3], routeId)
        seq = tmpIn[4]
        if int(seq)==1:
            prevNodeId = nodeId
            prevNodeTime = float(tmpIn[1])
        if int(seq)>1:
            
            dist = haversine(nodeSet[prevNodeId[0]].long, nodeSet[prevNodeId[0]].lat, nodeSet[nodeId[0]].long, nodeSet[nodeId[0]].lat)
            time = float(tmpIn[1])
            #print([lineSet[k].lineId for k in lineSet if lineSet[k].lineName == tripId])
            #Id = [lineSet[k].lineId for k in lineSet if tripId in lineSet[k].trips]
            if len([k for k in lineSet if tripId in lineSet[k].trips]) == 1:
                #Id = Id[0]
                linkId = (prevNodeId, nodeId)
            else:
                print("Multiple tripsIds are in the given route")
                print([k for k in lineSet if tripId in lineSet[k].trips])
                
            if prevNodeId not in nodeSet:
                nodeSet[prevNodeId] = Node([nodeSet[prevNodeId[0]].lat, nodeSet[prevNodeId[0]].long, prevNodeId[1]])
                
            if nodeId not in nodeSet:
                nodeSet[nodeId] = Node([nodeSet[nodeId[0]].lat, nodeSet[nodeId[0]].long, nodeId[1]])
                
            if linkId not in linkSet:
                linkSet[linkId] = Link([prevNodeId, nodeId, dist, (time - prevNodeTime)/60], "transit")
                linkSet[linkId].lineId = routeId
                if linkId not in nodeSet[prevNodeId].outLinks:
                    nodeSet[prevNodeId].outLinks.append(linkId)
                if linkId not in nodeSet[nodeId].inLinks:
                    nodeSet[nodeId].inLinks.append(linkId)
            linkId = (nodeId, prevNodeId)       
            if linkId not in linkSet:
                linkSet[linkId] = Link([nodeId, prevNodeId, dist, (time - prevNodeTime)/60], "transit")
                linkSet[linkId].lineId = routeId
                if linkId not in nodeSet[prevNodeId].inLinks:
                    nodeSet[prevNodeId].inLinks.append(linkId)
                if linkId not in nodeSet[nodeId].outLinks:
                    nodeSet[nodeId].outLinks.append(linkId) 
            
            prevNodeId = nodeId # Changing the previous node to current node (This will take care of the stop seq = 1 also)
            prevNodeTime = time
        else:
            continue
    inFile.close()

def readTransferLinks():
    # Reading transit transfer links
    transitLinks = [(k[0][0], k[1][0]) for k in linkSet if linkSet[k].type == 'transit']
    inFile = open(loc+"ft_input_transfers.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        if (tmpIn[0], tmpIn[1]) not in transitLinks and (tmpIn[1], tmpIn[0]) not in transitLinks:    
            linkId = (tmpIn[0], tmpIn[1])
            if linkId not in linkSet:
                linkSet[linkId] = Link(tmpIn, "transitTransfer")
                if linkId not in nodeSet[tmpIn[0]].outLinks:
                    nodeSet[tmpIn[0]].outLinks.append(linkId)
                if linkId not in nodeSet[tmpIn[1]].inLinks:
                    nodeSet[tmpIn[1]].inLinks.append(linkId)
                

    inFile.close()
    
def createBoardAlightLinks():
    transitNodes = list({k for k in nodeSet if nodeSet[k].type not in ['road', 'stop']})
    for l in transitNodes:
        linkId = (l[0], l)
        if linkId not in linkSet:
            linkSet[linkId] = Link([l[0], l, 0, 0], "boarding")
            if linkId not in nodeSet[l[0]].outLinks:
                    nodeSet[l[0]].outLinks.append(linkId)
            if linkId not in nodeSet[l].inLinks:
                    nodeSet[l].inLinks.append(linkId)
        linkId = (l, l[0])
        if linkId not in linkSet:
            linkSet[linkId] = Link([l, l[0], 0, 0], "alighting")
            if linkId not in nodeSet[l].outLinks:
                    nodeSet[l].outLinks.append(linkId)
            if linkId not in nodeSet[l[0]].inLinks:
                    nodeSet[l[0]].inLinks.append(linkId)
        
def readDemand():
    '''
    Reads passenger origin and destination
    '''
    inFile = open(loc + "demand.dat")
    tmpIn = inFile.readline().strip().split("\t")
    Id = 1
    for x in inFile:
        tmpIn = x.strip().split("\t")
        passengerSet[Id] = Passenger(tmpIn)
        Id = Id +1
        pairId = (tmpIn[0], tmpIn[1])
        if tmpIn[1] not in zoneSet[tmpIn[0]].dest:
            zoneSet[tmpIn[0]].dest.append(tmpIn[1])
            
        if pairId not in tripSet:
            tripSet[pairId] = Demand(tmpIn)
        else:
            print("O-D pair already present", tmpIn[0], tmpIn[1])
    inFile.close()    
    
def linkSetLines():
    '''
        Desfines the set of links associated with a particular line Id
    '''
    for l in linkSet:
        if linkSet[l].lineId != "":
            if l not in lineSet[linkSet[l].lineId].links:
                lineSet[linkSet[l].lineId].links.append(l)  
    
 ################################################################################################

def plotNetwork():
    roadNodes = [k for k in nodeSet if nodeSet[k].type == 'road']
    transitNodes = [k for k in nodeSet if nodeSet[k].type == 'transit']    
    roadLinks = [k for k in linkSet if linkSet[k].type == 'road']
    transitLinks = [k for k in linkSet if linkSet[k].type == 'transit']
    waitingLinks = [k for k in linkSet if linkSet[k].type in ['transitAccess', 'transitEgress', 'transitTransfer']]
    
    
    nodeTypes = ['Road nodes', 'Transit nodes']
    nodeColors = ['black', 'grey']
    nodeSizes = [500, 150]
    nodeTypeDict = dict(zip(nodeTypes, [roadNodes, transitNodes]))
    nodeColorDict = dict(zip(nodeTypes, nodeColors))
    nodeSizeDict = dict(zip(nodeTypes, nodeSizes))
    nodePos  = dict(zip(nodeSet,[(nodeSet[n].long,nodeSet[n].lat)
                                            for n in nodeSet]))
    
    
    
    linkTypes = ['Road links', 'Transit links', 'Waiting links']
    linkColors = ['black', 'grey', 'black']
    linkStyles = ['dashed', 'dashed', 'solid']
    linkTypeDict = dict(zip(linkTypes, [roadLinks, transitLinks, waitingLinks]))
    linkColorDict = dict(zip(linkTypes, linkColors))    
    linkSylesDict = dict(zip(linkTypes, linkStyles))



    g = nx.DiGraph()
    g.add_nodes_from(nodeSet)
    g.add_edges_from(roadLinks)

    fig, ax = plt.subplots(1, figsize=(40, 25))
    plt.axis('off')
    # iterate each nodetype, changing colors and labels of the nodes
    for nt in nodeTypes:
        # choose nodes and color for each iteration
        nlist = nodeTypeDict[nt]
        ncolor = nodeColorDict[nt]
        # draw the graph
        nx.draw_networkx_nodes(g,
                               pos=nodePos,
                               nodelist=nlist,
                               ax=ax,
                               node_color=ncolor,
                               label=nt, node_size = nodeSizeDict[nt])

    for k in linkTypes:
        nx.draw_networkx_edges(g, nodePos, edgelist=linkTypeDict[k], ax=ax, arrows=False, label=k, width = 2, edge_color = linkColorDict[k], style = linkSylesDict[k])
    #ed =nx.draw_networkx_edges(g, nodePos, edgelist=roadLinks, arrows=False, width = 2, edge_color = 'b')
    #ed =nx.draw_networkx_edges(g, nodePos, edgelist=waitingLinks, style = 'dashed', arrows=False, width = 2, edge_color = 'r')
    #ed = nx.draw_networkx_edges(g, nodePos, edgelist=roadLinks, edge_color=weights, arrows=False, width = 3, edge_cmap=plt.cm.YlOrRd)
    #ed =nx.draw_networkx_edges(g, nodePos, edgelist=transitLinks, edge_color=weights, arrows=False, width = 3, edge_cmap=plt.cm.YlOrRd)
    #nx.draw_networkx_edges(g, nodePos, edgelist=waitingLinks, arrows=False, width = 2, edge_color = 'r')
    #plt.colorbar(ed, ax=ax)
    #ed.set_clim(vmin=0, vmax=10)
    
    ax.legend(scatterpoints=1, loc='upper center', bbox_to_anchor=(0.5, 0.35, 0.5, 0.5), fontsize =  'xx-large')
    plt.show()
    
    

def analyzeTimeDistr():
    '''
        Analyzing travel time of different types of links so that the results are not biased
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    import random
    types = list({linkSet[l].type for l in linkSet})
    n_bins = 5


    for t in types:
        plt.hist([linkSet[k].time for k in linkSet if linkSet[k].type == t], n_bins,  histtype='bar', stacked=True, label=t)
    plt.ylabel('Frequency')
    plt.xlabel('Time (min)')
    plt.legend(loc="upper right")
    plt.show()
################################################################################################
def transitAssignment():
    '''
    This solves the transit assignment problem
    :return:
    '''
    
    m = Model()
    v = {(a, k): m.addVar(vtype=GRB.CONTINUOUS, obj =1.0, lb = 0.0, name='_'.join(['v', str(a), k])) for a in linkSet for k in zoneSet}
    Wt = {(i, k): m.addVar(vtype=GRB.CONTINUOUS, obj =1.0,  lb = -GRB.INFINITY , ub = GRB.INFINITY, name= '_'.join(['wt', i, k])) for i in stops for k in zoneSet}
    Wr = {(i, k): m.addVar(vtype=GRB.CONTINUOUS, obj =1.0,  lb = -GRB.INFINITY , ub = GRB.INFINITY, name= '_'.join(['wr', i, k])) for i in zoneSet for k in zoneSet}

    m.update()
    consConstr ={};rWaitConstr = {};tWaitConstr ={}
    for k in zoneSet:
        for i in nodeSet:
            tmp = sum([v[j, k] for j in nodeSet[i].outLinks]) - sum([v[j, k] for j in nodeSet[i].inLinks])
            if i == k:
                dem = sum([tripSet[d].demand for d in tripSet if d[1] == k])
                consConstr[i, k] = m.addConstr(tmp == -dem)
            elif (i, k) in tripSet:
                consConstr[i, k] = m.addConstr(tmp == tripSet[i, k].demand)
            else:
                consConstr[i, k] = m.addConstr(tmp == 0)

        for i in stops:
            for l in nodeSet[i].outLinks:
                if linkSet[l].type == 'boarding':
                    m.addConstr(v[l, k] <= (12/60.0)*Wt[i, k])

        for i in zoneSet: 
            for a in nodeSet[i].outLinks:
                if linkSet[a].type == 'road':
                    m.addConstr(v[a, k] <= A * 50*Wr[i, k])

    ivt = sum([v[l, k] * linkSet[l].time for l in linkSet for k in zoneSet])
    twt = sum([Wt[k] for k in Wt])
    rwt = sum([Wr[k] for k in Wr])
    obj = ivt +rwt + twt 
    m.setObjective(obj, sense=GRB.MINIMIZE)
    m.update()
    m.Params.OutputFlag = 0
    m.Params.DualReductions = 0
    m.optimize()
    print("Solving model took ", time.time() - start, " sec", m.status)
    print(m.objVal)
    print(ivt.getValue())
    print(twt.getValue())
    print(rwt.getValue())
    

def solveGurobiModel():
    m = Model()
    x = {l: m.addVar(vtype=GRB.BINARY, obj=1.0, lb=0.0, name=str(l)) for l in lineSet}
    y = {(l, f): m.addVar(vtype=GRB.BINARY, obj=1.0, lb=0.0, name='_'.join([l, str(f)])) for f in freqSet for l in lineSet}
    N = {(i, n): m.addVar(vtype=GRB.BINARY, obj=1.0, lb=0.0, name='_'.join([i, str(n)])) for i in zoneSet for n in fleetSet}
    

    v = {(a, k): m.addVar(vtype=GRB.CONTINUOUS, obj =1.0, lb = 0.0, name='_'.join(['v', str(a), k])) for a in linkSet for k in zoneSet}
    Wt = {(i, k): m.addVar(vtype=GRB.CONTINUOUS, obj =1.0,  lb = -GRB.INFINITY , ub = GRB.INFINITY, name= '_'.join(['wt', i, k])) for i in stops for k in zoneSet}
    Wr = {(i, k): m.addVar(vtype=GRB.CONTINUOUS, obj =1.0,  lb = -GRB.INFINITY , ub = GRB.INFINITY, name= '_'.join(['wr', i, k])) for i in zoneSet for k in zoneSet}
    t = {}; omega ={}
    for k in zoneSet:
        for i in stops:
            for a in nodeSet[i].outLinks:
                if linkSet[a].type == 'boarding':
                    for f in freqSet:
                        t[f, a, i, k] = m.addVar(vtype=GRB.CONTINUOUS, lb = 0.0, name= '_'.join(['t', str(f), str(a), i, k]))
        
        for i in zoneSet:
            for a in nodeSet[i].outLinks:
                if linkSet[a].type == 'road':
                    for n in fleetSet:
                        omega[n,a, i, k] = m.addVar(vtype=GRB.CONTINUOUS, lb = 0.0, name= '_'.join(['omega', str(n), str(a), i, k]))
    m.update()    
  
    
    tempVeh = 0
    for l in lineSet:
        tempTime = sum([linkSet[a].time for a in lineSet[l].links])
        tempVeh +=  tempTime * (sum([(f/60.0) * y[l, f] for f in freqSet]))
        m.addConstr(sum([y[l, f] for f in freqSet]) == x[l])
    m.addConstr(tempVeh <= maxBusFleet)
    
    tempVeh = 0
    for i in zoneSet:
        tempVeh += sum([n * N[i, n] for n in fleetSet])
        m.addConstr(sum([N[i, n] for n in fleetSet]) == 1)
    m.addConstr(tempVeh <= maxAVfleet)
       
    for k in zoneSet:
        dem = sum([tripSet[d].demand for d in tripSet if d[1] == k])
        for i in nodeSet:
            tmp = sum([v[j, k] for j in nodeSet[i].outLinks]) - sum([v[j, k] for j in nodeSet[i].inLinks])
            if i == k:
                dem = sum([tripSet[d].demand for d in tripSet if d[1] == k])
                m.addConstr(tmp == -dem)
            elif (i, k) in tripSet:
                m.addConstr(tmp == tripSet[i, k].demand)
            else:
                m.addConstr(tmp == 0)

        for i in stops:
            for a in nodeSet[i].outLinks:
                if linkSet[a].type == 'boarding':
                    l = a[1][1]
                    m.addConstr(v[a, k] <=   sum([(f/60.0) * t[f, a, i, k] for f in freqSet]))
                    for f in freqSet:
                        m.addConstr(Wt[i, k] - t[f, a, i, k] <= BigM*(1-y[l, f]), name = 'c_' + str(f)+ str(a)+str(i) + str(k))
                        m.addConstr(t[f, a, i, k] <= BigM*y[l, f], name = 'm_' + str(f)+ str(a)+str(i) + str(k))
                        m.addConstr(Wt[i, k] >= t[f, a, i, k])
                        

                        
                
        for i in zoneSet: 
            for a in nodeSet[i].outLinks:
                if linkSet[a].type == 'road':
                    m.addConstr(v[a, k] <= A * sum([n * omega[n, a, i, k] for n in fleetSet]))
                    for n in fleetSet:
                        bigM = dem * (1/ (A ))
                        m.addConstr(Wr[i, k] - omega[n, a, i, k] <= BigM*(1 - N[i, n]))
                        m.addConstr(omega[n, a, i, k] <= BigM*N[i, n])
                        m.addConstr(Wr[i, k] >= omega[n, a, i, k])

                    
    m.update()
    ivt = sum([v[l, k] * linkSet[l].time for l in linkSet for k in zoneSet])
    twt = sum([Wt[k] for k in Wt])
    rwt = sum([Wr[k] for k in Wr])
    obj = ivt + rwt + twt
    m.setObjective(obj, sense=GRB.MINIMIZE)
    '''
    for k in zoneSet:
        for i in stops:
            for a in nodeSet[i].outLinks:
                if linkSet[a].type == 'boarding':
                    l = a[1][1]
                    for f in freqSet:
                        m.setObjective(t[f, a, i, k], sense=GRB.MAXIMIZE)
                        m.Params.OutputFlag = 0
                        m.update()
                        m.optimize()
                        if m.status == 2:                            
                            print(m.objVal)
                        else:
                            print('Not possible', m.status)
    '''
                        
    m.update()
    m.Params.OutputFlag = 1
    m.Params.DualReductions  = 0
    start = time.time()
    m.optimize()
    if m.status == 2:            
        print("Solving model took ", time.time() - start, " sec")
        print(m.status)
        print(m.objVal)
        print(ivt.getValue())
        print(twt.getValue())
        print(rwt.getValue())
        print([k for k in x if x[k].x != 0])
        
        #roadShare =  sum([v[l].x for l in v if linkSet[l[0]].type == 'road'])
        #transitShare = sum([linkFlow[l].x for l in linkFlow if linkSet[l[0]].type == 'access'])
        #multiModShare = sum([linkFlow[l].x for l in linkFlow if linkSet[l[0]].type == 'modeTransAccess'])
            
    else:
        print("Extreme ray encountered!", m.status)
        m.computeIIS()
        m.write("Infeasible_model.ilp")
        
def setupMasterProblemModel(type = 'classic', verbose = 0):
    '''
        Sets up the gurobi model
        type = can be classic or multiple 
        verbose = 0: no verbose, 1: print the iterations of Gurobi, 2: print status and objective values,

    '''
    m = Model()  
    x = {l: m.addVar(vtype=GRB.BINARY, obj=1.0, lb=0.0, name=str(l)) for l in lineSet}
    y = {(l, f): m.addVar(vtype=GRB.BINARY, obj=1.0, lb=0.0, name='_'.join([l, str(f)])) for f in freqSet for l in lineSet}
    N = {(i, n): m.addVar(vtype=GRB.BINARY, obj=1.0, lb=0.0, name='_'.join([i, str(n)])) for i in zoneSet for n in fleetSet}
    eta = {k: m.addVar(vtype = GRB.CONTINUOUS, lb = -GRB.INFINITY , ub =GRB.INFINITY, name = 'eta_'+ str(k)) for k in zoneSet}
    
    tempVeh = 0
    for l in lineSet:
        tempTime = sum([linkSet[a].time for a in lineSet[l].links])
        tempVeh +=  tempTime * (sum([(f/60.0) * y[l, f] for f in freqSet]))
        m.addConstr(sum([y[l, f] for f in freqSet]) == x[l])
    m.addConstr(tempVeh <= maxBusFleet)
    
    tempVeh = 0
    for i in zoneSet:
        tempVeh += sum([n * N[i, n] for n in fleetSet])
        m.addConstr(sum([N[i, n] for n in fleetSet]) == 1)
    m.addConstr(tempVeh <= maxAVfleet)
    

        
    m.setObjective(sum([eta[k] for k in eta]), sense=GRB.MINIMIZE)
    m.update()
    m.Params.OutputFlag = verbose   
    if type == 'multiple':        
        m.Params.lazyConstraints = 1
        # Limit how many solutions to collect
        m.setParam(GRB.Param.PoolSolutions, 2)
        # Limit the search space by setting a gap for the worst possible solution
        # that will be accepted
        m.setParam(GRB.Param.PoolGap, 0.20)
        # do a systematic search for the k-best solutions
        m.setParam(GRB.Param.PoolSearchMode, 1)       
        #m.Params.lazyConstraints = 1
    m.update()
    return(m)


                
def BendersSubProblem(x, y, N, m, verbose = 0, type='classic'):
    '''
        Solves the Benders subproblems and applies cuts based on the type of the algorithm
        verbose = 0: no verbose, 1: print status and objective values, 2: print the iterations of Gurobi
        type = can be classic, aggregated, multiple cuts, 

        
    '''
    
    m1 = Model()
    v = {(a, k): m1.addVar(vtype=GRB.CONTINUOUS, obj =1.0, lb = 0.0, name='_'.join(['v', str(a), k])) for a in linkSet for k in zoneSet}
    Wt = {(i, k): m1.addVar(vtype=GRB.CONTINUOUS, obj =1.0,  lb = -GRB.INFINITY , ub = GRB.INFINITY, name= '_'.join(['wt', i, k])) for i in stops for k in zoneSet}
    Wr = {(i, k): m1.addVar(vtype=GRB.CONTINUOUS, obj =1.0,  lb = -GRB.INFINITY , ub = GRB.INFINITY, name= '_'.join(['wr', i, k])) for i in zoneSet for k in zoneSet}
    t = {}; omega ={}
    for k in zoneSet:
        for i in stops:
            for a in nodeSet[i].outLinks:
                if linkSet[a].type == 'boarding':
                    for f in freqSet:
                        t[f, a, i, k] = m1.addVar(vtype=GRB.CONTINUOUS, lb = 0.0, name= '_'.join(['t', str(f), str(a), i, k]))
        for i in zoneSet:
            for a in nodeSet[i].outLinks:
                if linkSet[a].type == 'road':
                    for n in fleetSet:
                        omega[n,a, i, k] = m1.addVar(vtype=GRB.CONTINUOUS, lb = 0.0, name= '_'.join(['omega', str(n), str(a), i, k]))
    m1.update()    
    
    consConstr = {}; constr1 = {}; constr2 = {}; constr3 = {}; constr4 = {}; constr5 = {}; constr6 = {}; constr7 = {}; constr8 = {}                 
    for k in zoneSet:
        dem = sum([tripSet[d].demand for d in tripSet if d[1] == k])

        for i in nodeSet:
            tmp = sum([v[j, k] for j in nodeSet[i].outLinks]) - sum([v[j, k] for j in nodeSet[i].inLinks])
            if i == k:
                dem = sum([tripSet[d].demand for d in tripSet if d[1] == k])
                consConstr[i, k] = m1.addConstr(tmp == -dem)
            elif (i, k) in tripSet:
                consConstr[i, k] = m1.addConstr(tmp == tripSet[i, k].demand)
            else:
                consConstr[i, k] = m1.addConstr(tmp == 0)
        
        for i in stops:
            for a in nodeSet[i].outLinks:
                if linkSet[a].type == 'boarding':
                    l = a[1][1]
                    
                    constr1[a, i, k] = m1.addConstr(v[a, k] <=   sum([round(f/60.0, 2) * t[f, a, i, k] for f in freqSet]))
                    for f in freqSet:
                        bigM = round(dem*60.0/f, 2)
                        '''
                        if dem*60/f > BigM:
                            print(True, dem*60/f, f)
                        '''
                        constr2[f, a, i, k] = m1.addConstr(Wt[i, k] - t[f, a, i, k] <= bigM*(1-y[l, f]))
                        constr3[f, a, i, k] = m1.addConstr(t[f, a, i, k] <= bigM*y[l, f])
                        constr4[f, a, i, k] = m1.addConstr(Wt[i, k] >= t[f, a, i, k])
        
        for i in zoneSet: 
            for a in nodeSet[i].outLinks:
                if linkSet[a].type == 'road':
                    constr5[i, a, k] = m1.addConstr(v[a, k] <= A * sum([n * omega[n, a, i, k] for n in fleetSet]))
                    for n in fleetSet:
                        bigM = dem * (1/ (A ))
                        constr6[n, a, i, k] = m1.addConstr(Wr[i, k] - omega[n, a, i, k] <= bigM*(1 - N[i, n]))
                        constr7[n, a, i, k] = m1.addConstr(omega[n, a, i, k] <= bigM*N[i, n])
                        constr8[n, a, i, k] = m1.addConstr(Wr[i, k] >= omega[n, a, i, k])
                                
                    
    m1.update()
    ivt = sum([v[l, k] * linkSet[l].time for l in linkSet for k in zoneSet])
    twt = sum([Wt[k] for k in Wt])
    rwt = sum([Wr[k] for k in Wr])
    obj = ivt + rwt + twt
    m1.setObjective(obj, sense=GRB.MINIMIZE)
    m1.update()
    m1.Params.OutputFlag = 0
    m1.Params.DualReductions  = 0
    start = time.time()
    m1.optimize()
    if m1.status == 2:   
        '''         
        print("Solving model took ", time.time() - start, " sec")
        print(m1.status)
        print(m1.objVal)
        print(ivt.getValue())
        print(twt.getValue())
        print(rwt.getValue())
        '''

        expr = 0.0;expr1= 0.0;expr2= 0.0;expr3= 0.0;expr4= 0.0
        for k in zoneSet:
            
            for i in nodeSet:
                tmp = sum([v[j, k] for j in nodeSet[i].outLinks]) - sum([v[j, k] for j in nodeSet[i].inLinks])
                if i == k:
                    dem = sum([tripSet[d].demand for d in tripSet if d[1] == k])
                    expr += consConstr[i, k].pi * (-dem)
                elif (i, k) in tripSet:
                    dem = tripSet[i, k].demand
                    expr += consConstr[i, k].pi * tripSet[i, k].demand
                else:
                    expr += 0
            
            for i in stops:
                for a in nodeSet[i].outLinks:
                    if linkSet[a].type == 'boarding':
                        l = a[1][1]                        
                        for f in freqSet:
                            bigM = round(dem*60.0/f, 2)
                            expr1 += constr2[f, a, i, k].pi * bigM*(1-m.getVarByName('_'.join([l, str(f)])))
                            expr2 += constr3[f, a, i, k].pi * bigM*m.getVarByName('_'.join([l, str(f)]))
            
            for i in zoneSet:     
                for a in nodeSet[i].outLinks:
                    if linkSet[a].type == 'road':
                        for n in fleetSet:
                            bigM = dem * (1 / (A ))
                            expr3 += constr6[n, a, i, k].pi * bigM * (1 - m.getVarByName('_'.join([i, str(n)])))
                            expr4 += constr7[n, a, i, k].pi * bigM * m.getVarByName('_'.join([i, str(n)]))
                            
        m.addConstr(sum([m.getVarByName('eta_'+str(k)) for k in zoneSet]) >= expr  + expr1 + expr2 + expr3 + expr4)   
        '''
        tmp = 0
        for k in zoneSet:
            for l in lineSet:
                for f in freqSet:
                    if y[l, f] > 0.4:
                        tmp +=  (1-m.getVarByName('_'.join([l, str(f)])))
                    else:
                        tmp +=  m.getVarByName('_'.join([l, str(f)]))                    
            
            for i in zoneSet:
                for n in fleetSet:            
                    if N[i, n] > 0.4:
                       tmp +=  (1 - m.getVarByName('_'.join([str(i), str(n)])))
                    else:
                        tmp +=  m.getVarByName('_'.join([str(i), str(n)]))
        for l in lineSet:
            if x[l] > 0.4:
                tmp +=  (1-m.getVarByName(str(l)))
            else:
                tmp +=  (m.getVarByName(str(l)))
                
        m.addConstr((m1.objVal-687555) * tmp + sum([m.getVarByName('eta_'+str(k)) for k in zoneSet]) >= m1.objVal)   

        m.update()
        '''
        return m, m1.objVal

    else:
        print("Extreme ray encountered!", m1.status)
        m1.computeIIS()
        m1.write("Infeasible_model.ilp")
                

def Benders(eps=1000, maxIt=1000, type= 'classic'):
    '''
    Implements benders decomposition
    eps = tolerance in the UB and LB
    maxIt = maximum iterations 
    type = can be classic, aggregated, multiple cuts, 
    
    '''
    UB = float("inf")
    LB = -float("inf")
    optCuts = []
    tol = float("inf")
    it = 0
    '''
    Computing initial values for the subproblem
    '''
    x0 = {l: 1 for l in lineSet}
    y0 = {(l, f): 1 for f in freqSet for l in lineSet}
    for k in y0:
        if k[1] != 12:
            y0[k] = 0
    N0 = {(i, n): 1 for i in zoneSet for n in fleetSet}
    for k in N0:
        if k[1] != 50:
            N0[k] = 0           
            
    m = setupMasterProblemModel(type=type)
    solutions = [(x0,y0,N0)]
    solutionsEverFound = [(x0,y0,N0)]
    
 
    while tol != 0 and it < maxIt:   
        for s in solutions:
            m, obj  = BendersSubProblem(s[0], s[1], s[2], m)
            
        UB = round(obj) #min(UB, obj)
        m.optimize()
        if m.status == 2:
            ob = m.objVal; x0 = {l:m.getVarByName(str(l)).x for l in lineSet}; y0 = {(l, f):m.getVarByName('_'.join([l, str(f)])).x  for f in freqSet for l in lineSet}; N0 =  {(i, n):m.getVarByName('_'.join([i, str(n)])).x  for i in zoneSet for n in fleetSet}
            LB = round(ob)
            solutions = []
            if (x0, y0, N0) not in solutionsEverFound :
                solutionsEverFound.append((x0, y0, N0))
                solutions.append((x0, y0, N0))
                
            if type == 'multiple':                
                nSolutions = m.SolCount            
                if (nSolutions >= 2):
                    for solNum in range(nSolutions - 1):                    
                        m.setParam(GRB.Param.SolutionNumber, solNum+1)
                        ob = m.objVal; x0 = {l:m.getVarByName(str(l)).x for l in lineSet}; y0 = {(l, f):m.getVarByName('_'.join([l, str(f)])).x  for f in freqSet for l in lineSet}; N0 =  {(i, n):m.getVarByName('_'.join([i, str(n)])).x  for i in zoneSet for n in fleetSet}
                        LB = max(LB, round(ob))
                        if (x0, y0, N0) not in solutionsEverFound:
                            solutionsEverFound.append((x0, y0, N0))
                            solutions.append((x0, y0, N0))
        else:
            print('master problem cannot be solved... \n Terminating ....')
            break
       
        tol = UB - LB
        it += 1
        print((it, LB, UB, tol,  [k for k in x0 if round(x0[k]) == 0]))
    return (x0, y0, N0)


            
            
            

    
################################################################################################
loc = 'Z:/Projects/NSF_SCC/Transit network design for FMLM in case of AVs/Transit FMLM AV/Scripts/InputFiles/Siuox Falls network/'
start = time.time()
alpha = [4, 2] # Fuel, transfer penalty
baseTaxiFare = 0.8 # in dollars
fuelCost = 0.21 # in dollars
transitFare = 2 # in dollars applies to only access and mode transfer links
VOT = 23 # in dollars per minute
freqSet =[12, 6]# [2, 3, 4, 6, 12] #  Buses per hour
#freqSet = [1/60, 1/30, 1/6, 1/4, 1/2]
#freqSet = [float(f) for f in freqSet]
fleetSet = [1, 50] # [1, 50, 100, 200, 500] # 
tBigM = 200
BigM = 100000
rWaitFac = 1000
A =  6.06/3600 #48.31 # 6.06 for actual alpha

maxBusFleet = 200
maxAVfleet = 3000


zoneSet = {}
nodeSet = {}
linkSet ={}
lineSet ={}
passengerSet={}
tripSet ={}
transitWaitNodeDict = defaultdict(list)
transitWaitingNodes = []
waitingLineFinder = {}





readNodes()
readLinks()
readLines()
readtrips()
readTransitLinks()
readTransferLinks()
linkSetLines()
createBoardAlightLinks()
readDemand()
print(len(nodeSet), "nodes in the network")
print(len(linkSet), "links in the network")
print("Reading network took ", round(time.time() - start), " seconds")
stops = list({k for k in nodeSet if nodeSet[k].type == 'stop'})

#############################################################################################
# Testing
start = time.time()
#solveGurobiModel()
Benders()
print('It took ', time.time() - start)




