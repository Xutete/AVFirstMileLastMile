# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 18:37:08 2020
This serves as the final code for the AV-TNDP design. 
I'll try to add more and more details in the script. 


@author: Pramesh Kumar
"""

import math, time
from gurobipy import *
from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

#########################################################################################################

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
            

########################################################################################################
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

#############################################################################################################
loc = 'S:/Projects/NSF_SCC/Transit network design for FMLM in case of AVs/Transit FMLM AV/Scripts/InputFiles/Siuox Falls network/'
start = time.time()
alpha = [4, 2] # Fuel, transfer penalty
baseTaxiFare = 0.8 # in dollars
fuelCost = 0.21 # in dollars
transitFare = 2 # in dollars applies to only access and mode transfer links
VOT = 23 # in dollars per minute
freqSet = [2, 3, 4, 6, 12] # [2,12]#  Buses per hour
#freqSet = [1/60, 1/30, 1/6, 1/4, 1/2]
#freqSet = [float(f) for f in freqSet]
fleetSet =  [1, 50, 100, 200, 500] # [50, 100] # 
tBigM = 100
BigM = 100000
rWaitFac = 1000
A =  6.06/3600 #48.31 # 6.06 for actual alpha

maxBusFleet = 100
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
'''
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
'''
