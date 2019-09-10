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
        self.lat = float(_tmpIn[1])
        self.long = float(_tmpIn[2])
        self.type = _timpIn[3]
        self.nodes = []



################################################################################################
def readNodes():
    # Zone file has "zoneId, Latitude, Longitude"
    inFile = open(inputDataLocation+"ft_input_zones.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        zoneId = tmpIn[0]

        zoneSet[zoneId] = Zone(tmpIn)
    inFile.close()
    print(len(zoneSet), "zones")

def readStops():
    inFile = open(inputDataLocation+"ft_input_stops.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        stopSet[tmpIn[0]] = StopNode(tmpIn)
    inFile.close()
    print(len(stopSet), "stops")

def readRoadNodes():
    inFile = open(inputDataLocation+"ft_input_roadNodes.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        stopSet[tmpIn[2]] = RoadNode(tmpIn)
    inFile.close()
    print(len(stopSet), "road nodes")


################################################################################################
zoneSet = {}
stopSet = {}
roadNodeSet ={}

readZones()
readStops()
readNodes()