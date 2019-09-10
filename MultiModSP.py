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
        self.nodes = []
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
    # Reading zones
    inFile = open(inputDataLocation+"ft_input_zones.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        zoneId = tmpIn[0]
        if zoneId not in nodeSet:
            nodeSet[zoneId] = Node([tmpIn[1], tmpIn[2], "zone"])
    inFile.close()

    # Reading transit stops
    inFile = open(inputDataLocation+"ft_input_stops.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        if tmpIn[0] not in nodeSet:
            nodeSet[tmpIn[0]] = Node([tmpIn[3], tmpIn[4], "transit"])
    inFile.close()

    # Reading road nodes
    inFile = open(inputDataLocation+"ft_input_roadNodes.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        if tmpIn[2] not in nodeSet:
            nodeSet[tmpIn[2]] = Node([tmpIn[0], tmpIn[1], "road"])
    inFile.close()
    print(len(nodeSet), "nodes in the network")

def readLinks():
    inFile = open(inputDataLocation+"ft_input_accessLinks.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        if (tmpIn[0], tmpIn[1]) not in linkSet:
            linkSet[tmpIn[0], tmpIn[1]] = Link([tmpIn[0], tmpIn[1], tmpIn[2], tmpIn[3], "access"])
        if (tmpIn[1], tmpIn[0]) not in linkSet:
            linkSet[tmpIn[0], tmpIn[1]] = Link([tmpIn[1], tmpIn[0], tmpIn[2], tmpIn[3], "egress"])
    inFile.close()

    print(len(linkSet), "links in the network")





################################################################################################
nodeSet = {}
linkSet ={}



readNodes()
readLinks()