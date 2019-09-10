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
        self.labels = (999999.0,999999.0) # time, cost
        self.preds = ("","")




class Link:
    def __init__(self, _tmpIn):
        self.fromNode = _tmpIn[0]
        self.toNode = _tmpIn[1]
        self.dist = _tmpIn[2] # in miles
        self.time = _tmpIn[3] # in minutes
        self.type = _tmpIn[4]
        self.passengers = []



################################################################################################
def readNodes():
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
    # Reading access links
    inFile = open(inputDataLocation+"ft_input_accessLinks.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        if (tmpIn[0], tmpIn[1]) not in linkSet:
            linkSet[tmpIn[0], (tmpIn[1], "transit")] = Link([tmpIn[0], tmpIn[1], tmpIn[2], tmpIn[3], "access"])
            nodeSet[tmpIn[0]].outLinks.append((tmpIn[1], "transit"))
            nodeSet[(tmpIn[1], "transit")].inLinks.append(tmpIn[0])
        if (tmpIn[1], tmpIn[0]) not in linkSet:
            nodeSet[(tmpIn[1], "transit")].outLinks.append()
            linkSet[tmpIn[0], tmpIn[1]] = Link([tmpIn[1], tmpIn[0], tmpIn[2], tmpIn[3], "egress"])
    inFile.close()

    # Reading
    print(len(linkSet), "links in the network")





################################################################################################
nodeSet = {}
linkSet ={}



readNodes()
readLinks()