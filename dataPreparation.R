##############################################
# This program creates input files for creating multimodal network
##############################################
# This is the query for the GTFS data. Please note that this is Green Line GTFS data for eastbound direction
library(RPostgreSQL)
drv <- dbDriver("PostgreSQL")
library(multiplex)
location <- "S:\\Projects\\NSF_SCC\\Transit network design for FMLM in case of AVs\\Transit FMLM AV\\Scripts\\InputFiles\\"
# Connecting to the server 
con <- dbConnect(drv, dbname = "gtfs",
                 host = "localhost", port = 9876,
                 user = "postgres", password = "postgres00")
#Querying the serviceID within your time period
serviceId <- dbGetQuery(con, "select service_id from gtfs_feb_mar_2016.calendar where start_date<'2016-03-08' and end_date>'2016-03-08'")
#Filtering for typical weekday service
serviceId <- serviceId[grepl('Weekday-01', serviceId$service_id) == TRUE & grepl('01-', serviceId$service_id) == FALSE, ]

# Querying and combining the gtfs data for different serviceIds (e.g. bus, light rail, sub, etc.)
gtfs <- data.frame()
for (sId in serviceId){
  gtfs <- rbind(gtfs, dbGetQuery(con, paste0("select start_date, end_date, routes.route_id, route_type, route_short_id, route_long_name, trips.trip_id, trips.direction_id, trip_headsign, arrival_time, departure_time, stops.stop_id, stop_name, stop_sequence, stop_desc, stop_lat, stop_lon from gtfs.gtfs_feb_mar_2016.calendar inner join gtfs.gtfs_feb_mar_2016.trips ON (calendar.service_id = trips.service_id and (calendar.service_id = ", "'", sId, "'", ")) join gtfs.gtfs_feb_mar_2016.stop_times ON (stop_times.trip_id = trips.trip_id) join gtfs.gtfs_feb_mar_2016.routes ON (routes.route_id = trips.route_id) join gtfs.gtfs_feb_mar_2016.stops ON (stops.stop_id = stop_times.stop_id)")))
}

#Changing the format of time for ft_input
gtfs$arrival_time <- gsub(":", "", as.character(gtfs$arrival_time))
gtfs$departure_time <- gsub(":", "", as.character(gtfs$departure_time))




# Extracting trips on routes
trips <- gtfs[gtfs$stop_sequence == 1, ]
trips <- trips[c("trip_id", "route_id", "route_type", "departure_time", "direction_id")]
trips$capacity <- ""
trips$shapeId <- ""
trips<- trips[, c("trip_id", "route_id", "route_type", "departure_time", "capacity", "shapeId", "direction_id")]


# Extracting the number of stops for each trip
uniqueRoutes <- unique(trips$route_id)
trips$numStops <- 0
for (i in uniqueRoutes){
  selectedTrips <- trips[trips$route_id == i, ]
  for(j in seq(from = 1, to = nrow(selectedTrips))){
    trips[trips$trip_id == selectedTrips$trip_id[j], ]$numStops<- nrow(gtfs[gtfs$trip_id == selectedTrips$trip_id[j], ])
  }
}

# Extracting trips to be extracted from gtfs

extrTrips <- c()
for (i in uniqueRoutes){
  p <- trips[trips$route_id == i, ]
  ind <- which.max(p$numStops)
  extrTrips <- c(extrTrips, p[ind, ]$trip_id)
}


p <- gtfs[gtfs$trip_id %in% extrTrips, ]


# Writing the ft_input_route file
routes <- p[c('route_id', 'route_short_id', 'route_long_name', 'route_type')]
routes <- unique(routes)
write.table(routes, paste0(location, "ft_input_routes.dat"), sep = "\t", row.names = FALSE, quote = FALSE)

# Writing the ft_input_stops.dat
stops <- p[c("stop_id", "stop_name", "stop_desc", "stop_lat", "stop_lon")]
stops <- unique(stops)
stops$capacity <- 100
write.table(stops, paste0(location, "ft_input_stops.dat"), sep = "\t", row.names = FALSE, quote = FALSE)


# Writing the ft_input_trips.dat
#trips <- p[p$stop_sequence == 1, ]
trips <- p[c("trip_id", "route_id", "route_type", "departure_time", "direction_id", "stop_id")]
trips$capacity <- ""
trips$shapeId <- ""
trips<- trips[, c("trip_id", "route_id", "route_type", "departure_time", "capacity", "shapeId", "direction_id", "stop_id")]
write.table(trips, paste0(location, "ft_input_trips.dat"), sep = "\t", row.names = FALSE, quote = FALSE)


# Writing the ft_input_stopTimes.dat
stopTimes <- p[c("trip_id", "arrival_time", "departure_time", "stop_id", "stop_sequence")]
stopTimes <- unique(stopTimes)
stopTimes <- stopTimes[order(stopTimes$trip_id, stopTimes$stop_sequence), ]
write.table(stopTimes, paste0(location, "ft_input_stopTimes.dat"), sep = "\t", row.names = FALSE, quote = FALSE)


# Writing the ft_input_zones.dat
location2 <- "S:\\Projects\\UPassClustering\\Data\\TAZ2010\\"
zones <- read.csv(paste0(location2, "TAZ2010WGSCentroid.csv"))
zones <- zones[c('TAZ', 'Y', 'X')]
colnames(zones) <- c('zoneId',	'Latitude',	'Longitude')
write.table(zones, paste0(location, "ft_input_zones.dat"), sep = "\t", row.names = FALSE, quote = FALSE)


# Writing the ft_input_accessLinks.dat
# This step requires the ft_input_stops.dat and ft_input_zones.dat
start_time <- Sys.time()
library(geosphere)
zones <- read.table(paste0(location, "ft_input_zones.dat"), header = TRUE)
stops <-read.table(paste0(location, "ft_input_stops.dat"), header = TRUE)
access_links <- data.frame()
colnames(access_links) <- c('TAZ',	'stopId', 'dist',	'time')
for (i in seq(from = 1, to = nrow(zones))){
  tmpStops <- stops
  tmpStops$dist <- 0.000621371*distHaversine(c(as.numeric(zones$Longitude[i]), as.numeric(zones$Latitude[i])), stops[c('stop_lon', 'stop_lat')])
  tmpStops <- tmpStops[tmpStops$dist <= 0.75, ]
  if (nrow(tmpStops) != 0){
    tmpStops$time <- tmpStops$dist*60/3
    tmpStops$TAZ <- zones$zoneId[i]
    tmpStops <- tmpStops[c('TAZ',	'stop_id', 'dist',	'time')]
    colnames(tmpStops) <- c('TAZ', 'stopId', 'dist',	'time')
    access_links <- rbind(access_links, tmpStops)
  }
}
write.table(access_links, paste0(location, "ft_input_accessLinks.dat"), sep = "\t", row.names = FALSE, quote = FALSE)
Sys.time() -start_time


# Writing the ft_input_transfers.dat
# This step requires the ft_input_stops.dat
library(geosphere)
stops <-read.table(paste0(location, "ft_input_stops.dat"), header = TRUE)
# Note that the second stop should not be the last stop condtion is not included. This will be taken care by the SBSP
createTransferLinks <- function(i, toStop){
  fromStop <- toStop[i, ]
  toStop$dist <- 0.000621371*distHaversine(c(as.numeric(fromStop$stop_lon), as.numeric(fromStop$stop_lat)), toStop[c('stop_lon', 'stop_lat')])
  toStop <- toStop[toStop$dist <= 0.25, ]
  toStop$fromStop <- fromStop$stop_id
  toStop$time <- toStop$dist*60/3
  toStop <- toStop[toStop$fromStop != toStop$stop_id, ]
  toStop <- toStop[c('fromStop',	'stop_id', 'dist',	'time')]
  toStop <- unique(toStop)
  colnames(toStop) <- c('fromStop', 'toStop', 'dist',	'time')
  toStop <- na.omit(toStop)
  return(toStop)
}

start_time <- Sys.time()
transfers_links <- do.call(rbind,lapply(1:nrow(stops),function(x) createTransferLinks(x, stops)))
Sys.time() -start_time
write.table(transfers_links, paste0(location, "ft_input_transfers.dat"), sep = "\t", row.names = FALSE, quote = FALSE)


# Preparing ft_drivingTime.dat
DrivingTime <- read.csv("S:\\Projects\\RideSharing Work\\Scripts\\travel_times.csv", header = FALSE)
colnames(DrivingTime) <- c("fromTAZ", "toTAZ", "travelTime")
write.table(DrivingTime, paste0(location, "ft_input_drivingTime.dat"), sep = "\t", row.names = FALSE, quote = FALSE)


# Writing the ft_input_demand.dat
# This step requires the deamnd data
location2 <- "S:\\Projects\\RideSharing Work\\Scripts\\Optimization Trials\\"
OD.trips <- read.csv(paste0(location2, "dem.csv"))
OD.trips$PDT <- round(OD.trips$dep_time/60)
for (i in seq(from = 1, to = nrow(OD.trips))){
  if (OD.trips$role[i] == "R"){
    OD.trips$PAT[i] = round(OD.trips$dep_time[i]/60) + 1.6*DrivingTime[DrivingTime$fromTAZ == OD.trips$ORIGZONE[i] & DrivingTime$toTAZ == OD.trips$DESTZONE[i], ]$travelTime
  }
  else{
    OD.trips$PAT[i] = round(OD.trips$dep_time[i]/60) + 1.6*DrivingTime[DrivingTime$fromTAZ == OD.trips$ORIGZONE[i] & DrivingTime$toTAZ == OD.trips$DESTZONE[i], ]$travelTime
  }
}

OD.trips$timePeriod <- "AM" 
OD.trips$passengerId<- rownames(OD.trips)
OD.trips$direction <- 1
OD.trips <- OD.trips[c("passengerId", "ORIGZONE", "DESTZONE", "TRIPMODE", "timePeriod", "direction", "PDT", "PAT", "role")]
colnames(OD.trips) <- c("passengerId", "OrigTAZ", "DestTAZ", "mode", "timePeriod", "direction", "PDT", "PAT", "role")



# Writing demand scenarios
location4 <- "S:\\Projects\\NSF_SCC\\Transit network design for FMLM in case of AVs\\Transit FMLM AV\\Scripts\\InputFiles\\Demand Scenarios\\"
dem <- OD.trips <- read.csv(paste0(location4, "FullDemand.csv"))
colnames(dem) <- c("OriginTAZ", "DestTAZ", "HHID")
# Writing ft_scenarios for 10 scenarios
for(i in seq(from = 1, to = 10)){
  sampleDemand <- dem[sample(nrow(dem), 10000), ]
  write.table(sampleDemand, paste0(location4, "Scenario", i, ".dat"), sep = "\t", row.names = FALSE, quote = FALSE)
}


# Writing the ft_input_roadNodes.dat
location5 <- "S:\\Projects\\Ridesharing Work\\Scripts\\TC network\\"
nodes <- read.csv(paste0(location5, "nodes.csv"))
nodes <- nodes[c("X", "Y", "N")]
write.table(nodes, paste0(location, "ft_input_roadNodes.dat"), sep = "\t", row.names = FALSE, quote = FALSE)
  

# Writing the ft_input_roadLinks.dat
library("readxl")
location3 <- "S:\\Projects\\Ridesharing Work\\OD files\\"
road <- read_excel(paste0(location3, "tc_network.xlsx"))
road <- road[c("tail node", "head node", "capacity (veh/h)", "fftt(min)", "length (miles)")]
colnames(road) <- c("tailNode", "headNode", "cap", "fft", "len")
write.table(road, paste0(location, "ft_input_roadLinks.dat"), sep = "\t", row.names = FALSE, quote = FALSE)

# Writing the ft_input_moadTransfer.dat
# This file requires ft_input_roadNodes.dat and ft_input_stops
start_time <- Sys.time()
library(geosphere)
nodes <- read.table(paste0(location, "ft_input_roadNodes.dat"), header = TRUE)
stops <-read.table(paste0(location, "ft_input_stops.dat"), header = TRUE)
modeTransferLink <- data.frame(matrix(0, ncol=4))
colnames(modeTransferLink) <- c('Node',	'stopId', 'dist',	'time')
modeTransferLink <- modeTransferLink[-1, ]
for (i in seq(from = 1, to = nrow(nodes))){
  tmpStops <- stops
  tmpStops$dist <- 0.000621371*distHaversine(c(as.numeric(nodes$X[i]), as.numeric(nodes$Y[i])), stops[c('stop_lon', 'stop_lat')])
  tmpStops <- tmpStops[tmpStops$dist <= 0.75, ]
  if (nrow(tmpStops) != 0){
    tmpStops$time <- tmpStops$dist*60/3
    tmpStops$Node <- nodes$N[i]
    tmpStops <- tmpStops[c('Node',	'stop_id', 'dist',	'time')]
    colnames(tmpStops) <- c('Node', 'stopId', 'dist',	'time')
    modeTransferLink <- rbind(modeTransferLink, tmpStops)
  }
}
write.table(modeTransferLink, paste0(location, "ft_input_moadTransfer.dat"), sep = "\t", row.names = FALSE, quote = FALSE)
Sys.time() -start_time




