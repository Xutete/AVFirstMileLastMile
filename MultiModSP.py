# -*- coding: utf-8 -*-
################################################################################################
### Multimodal shortest path algorithm ###
################################################################################################
'''
Code primarily written by Pramesh Kumar
Contact: kumar372@umn.edu
-------------------------------------------------------
'''
################################################################################################
################################################################################################
################################################################################################

import math, time, heapq

inputDataLocation = "S:/Projects/NSF_SCC\Transit network design for FMLM in case of AVs/Transit FMLM AV/Scripts/InputFiles/"

################################################################################################
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
    def __init__(self, _tmpIn):
        self.fromNode = _tmpIn[0]
        self.toNode = _tmpIn[1]
        self.dist = float(_tmpIn[2]) # in miles
        self.time = float(_tmpIn[3]) # in minutes
        self.type = _tmpIn[4]
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
        self.freq = 10 # Buses per hour

class Passenger:
    def __init__(self, _tmpIn):
        self.origin = _tmpIn[0]
        self.dest = _tmpIn[1]
        self.path = []




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
    # Reading transit stops
    inFile = open(inputDataLocation+"ft_input_stops.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        nodeId =  (tmpIn[0], "transit") # A word transit is added because there are a lot of similar roadId and stopId
        if nodeId not in nodeSet:
            nodeSet[nodeId] = Node([tmpIn[3], tmpIn[4], "transit"])
        else:
            print(tmpIn[0], " stop already present as ", nodeSet[tmpIn[0]].type)
    inFile.close()

    # Reading road nodes
    inFile = open(inputDataLocation+"ft_input_roadNodes.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        nodeId =  tmpIn[2]
        if nodeId not in nodeSet:
            nodeSet[nodeId] = Node([tmpIn[0], tmpIn[1], "road"])
        else:
            print(tmpIn[2], " roadNode already present as ", nodeSet[tmpIn[2]].type)
    inFile.close()
    print(len(nodeSet), "nodes in the network")

    # Reading zones
    inFile = open(inputDataLocation+"ft_input_zones.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        zoneId = tmpIn[0]
        if zoneId not in nodeSet:
            nodeSet[zoneId] = Node([tmpIn[1], tmpIn[2], "zone"])
    inFile.close()



def readLinks():
    """
    Read access, mode transfer, road, and
    transit transfer links
    """
    # Reading access links
    inFile = open(inputDataLocation+"ft_input_accessLinks.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        if (tmpIn[0], (tmpIn[1], "transit")) not in linkSet:
            linkSet[tmpIn[0], (tmpIn[1], "transit")] = Link([tmpIn[0], tmpIn[1], tmpIn[2], tmpIn[3], "access"])
            nodeSet[tmpIn[0]].outLinks.append((tmpIn[1], "transit"))
            nodeSet[(tmpIn[1], "transit")].inLinks.append(tmpIn[0])
        if ((tmpIn[1], "transit"), tmpIn[0]) not in linkSet:
            linkSet[(tmpIn[1], "transit"), tmpIn[0]] = Link([tmpIn[1], tmpIn[0], tmpIn[2], tmpIn[3], "egress"])
            nodeSet[(tmpIn[1], "transit")].outLinks.append(tmpIn[0])
            nodeSet[tmpIn[0]].inLinks.append((tmpIn[1], "transit"))
    inFile.close()

    # Reading mode transfer links
    inFile = open(inputDataLocation+"ft_input_moadTransfer.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")

        if (tmpIn[0], tmpIn[1]) not in linkSet:
            linkSet[tmpIn[0], (tmpIn[1], "transit")] = Link([tmpIn[0], tmpIn[1], tmpIn[2], tmpIn[3], "modeTransAccess"])
            nodeSet[tmpIn[0]].outLinks.append((tmpIn[1], "transit"))
            nodeSet[(tmpIn[1], "transit")].inLinks.append(tmpIn[0])
        if ((tmpIn[1], "transit"), tmpIn[0]) not in linkSet:
            linkSet[(tmpIn[1], "transit"), tmpIn[0]] = Link([tmpIn[1], tmpIn[0], tmpIn[2], tmpIn[3], "modeTransEgress"])
            nodeSet[(tmpIn[1], "transit")].outLinks.append(tmpIn[0])
            nodeSet[tmpIn[0]].inLinks.append((tmpIn[1], "transit"))
    inFile.close()

    # Reading road links
    inFile = open(inputDataLocation+"ft_input_roadLinks.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        if (tmpIn[0], tmpIn[1]) not in linkSet:
            linkSet[tmpIn[0], tmpIn[1]] = Link([tmpIn[0], tmpIn[1], tmpIn[2], tmpIn[3], "road"])
            if tmpIn[0] not in nodeSet:
                nodeSet[tmpIn[0]] = Node([0, 0, "road"])
            if tmpIn[1] not in nodeSet:
                nodeSet[tmpIn[1]] = Node([0, 0, "road"])
            nodeSet[tmpIn[0]].outLinks.append(tmpIn[1])
            nodeSet[tmpIn[1]].inLinks.append(tmpIn[0])

    inFile.close()

    # Reading transit transfer links
    inFile = open(inputDataLocation+"ft_input_transfers.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        tailNode = (tmpIn[0], "transit")
        headNode = (tmpIn[1], "transit")
        if (tailNode, headNode) not in linkSet:
            linkSet[tailNode, headNode] = Link([tmpIn[0], tmpIn[1], tmpIn[2], tmpIn[3], "transitTransfer"])
            if tailNode not in nodeSet:
                nodeSet[tailNode] = Node([0, 0, "transit"])
            if headNode not in nodeSet:
                nodeSet[headNode] = Node([0, 0, "transit"])
            nodeSet[tailNode].outLinks.append(headNode)
            nodeSet[headNode].inLinks.append(tailNode)
    inFile.close()


def readLines():
    # Reading transit stops
    inFile = open(inputDataLocation + "ft_input_routes.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        lineSet[tmpIn[0]] = Line([tmpIn[0], tmpIn[2]])
    inFile.close()

def readtrips():
    # Reading transit stops
    inFile = open(inputDataLocation + "ft_input_trips.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        lineSet[tmpIn[1]].trips.append(tmpIn[0])
        lineSet[tmpIn[1]].stops.append(tmpIn[7])
    inFile.close()




def readTransitLinks():
    inFile = open(inputDataLocation+"ft_input_stopTimes.dat")
    tmpIn = inFile.readline().strip().split("\t")
    prevNodeId = ""
    for x in inFile:
        tmpIn = x.strip().split("\t")
        tripId = tmpIn[0]
        nodeId = (tmpIn[3], "transit")
        seq = tmpIn[4]
        if int(seq)==1:
            prevNodeId = nodeId
            prevNodeTime = float(tmpIn[1])
        if int(seq)>1:
            linkId = (prevNodeId, nodeId)
            dist = haversine(nodeSet[prevNodeId].long, nodeSet[prevNodeId].lat, nodeSet[nodeId].long, nodeSet[nodeId].lat)
            time = float(tmpIn[1])
            #print([lineSet[k].lineId for k in lineSet if lineSet[k].lineName == tripId])
            Id = [lineSet[k].lineId for k in lineSet if tripId in lineSet[k].trips]
            if len(Id) == 1:
                Id = Id[0]
            else:
                print("Multiple tripsIds are in the given route")
            linkSet[linkId] = Link([prevNodeId, nodeId, dist, (time - prevNodeTime), tripId])
            linkSet[linkId].lineId = Id
            nodeSet[prevNodeId].outLinks.append(nodeId)
            nodeSet[nodeId].inLinks.append(prevNodeId)
            prevNodeId = nodeId # Changing the previous node to current node (This will take care of the stop seq = 1 also)
            prevNodeTime = time
    inFile.close()
    print(len(linkSet), "links in the network")

def creatingUniqueSet():
    '''
    Create unique set of neighbours for each node in the network
    '''
    for i in nodeSet:
        nodeSet[i].inLinks = list(set(nodeSet[i].inLinks))
        nodeSet[i].outLinks = list(set(nodeSet[i].outLinks))

    for j in lineSet:
        lineSet[j].stops = list(set(lineSet[j].stops))
        lineSet[j].trips = list(set(lineSet[j].trips))



def readDemand():
    '''
    Reads passenger origin and destination
    '''
    inFile = open(inputDataLocation + "Demand Scenarios/Scenario1.dat")
    tmpIn = inFile.readline().strip().split("\t")
    Id = 1
    for x in inFile:
        tmpIn = x.strip().split("\t")
        passengerSet[Id] = Passenger(tmpIn)
        Id = Id +1
    inFile.close()





################################################################################################

def calcWaitFuelTime():
    '''
    Calculate the wait time of the links if they are access,
    transfer, or mode transfer links

    Calcualte the fuel cost for road links
    '''

    for l in linkSet:
        if linkSet[l].type in ['access', 'modeTransAccess', 'transitTransfer']:
            lineFreq = [lineSet[k].freq for k in lineSet if l[1][0] in lineSet[k].stops]
            linkSet[l].wait = 60.0/(sum(lineFreq))
        if linkSet[l].type == "road":
            linkSet[l].fuelCost = linkSet[l].dist * alpha[0] * fuelCost / 60




def shortestPath(origin):
    '''
    Calculate the shoetest path from a node to all other nodes
    '''
    label = {}
    pred ={}
    for n in nodeSet:
        nodeSet[n].label = float("inf")
        nodeSet[n].pred = ""
    nodeSet[origin].label = 0.0
    nodeSet[origin].pred = "NA"
    SE = [(0.0, origin)]
    while SE:
        i = heapq.heappop(SE)[1]
        #i = min(SE, key=label.get)
        #SE.remove(i)
        for j in nodeSet[i].outLinks:
            tmpCost = linkSet[i, j].time + linkSet[i, j].wait
            if linkSet[i, j].type == 'road':
                tmpCost = linkSet[i, j].dist*alpha[0]*fuelCost/60
            if nodeSet[j].label > nodeSet[i].label + tmpCost:
                nodeSet[j].label = nodeSet[i].label + tmpCost
                nodeSet[j].pred = i
                heapq.heappush(SE, (nodeSet[j].label , j))
                #SE.append(j)



def assignPassengers():
    for p in passengerSet:
        start = time.time()
        shortestPath(passengerSet[p].origin)
        dest = passengerSet[p].dest
        path=[dest]
        while nodeSet[dest].pred != "NA":
            path.append(linkSet[dest, nodeSet[dest].pred].type)
            dest = nodeSet[dest].pred
        print("Path of the passenger is ", path)
        print("One shortest path took ", time.time() - start, " seconds")
        print(nodeSet[passengerSet[p].dest].label)


################################################################################################
start = time.time()
nodeSet = {}
linkSet ={}
lineSet ={}
passengerSet={}
alpha = [15] #Fuel
fuelCost = 0.15


#if 'key' in myDict: del myDict['key']
readNodes()
readLinks()
readLines()
readtrips()
readTransitLinks()
creatingUniqueSet()
readDemand()
print("Creating the network and reading the demand took ", time.time() - start, " seconds")
calcWaitFuelTime()
start = time.time()
#assignPassengers()
print("assigning Passengers on shortest path is ", time.time() - start, " seconds")


