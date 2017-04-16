# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# THESIS: FIRST SEGMENTATION LEVEL - PART 2 ------------------------------------------------------------------------
#
# Raphael Knevels
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# DESCRIPTION:
# In this script the object-oriented image analysis for segmentation level 1
# is performed. The objective of this levels is the extraction of potential
# scarp objects out of the fine-scale segmentation.
#
# This script builds on the computated parameters in L1_part1.R.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



# CONTENT -----------------------------------------------------------------
# 1 PACKAGES, FUNCTIONS & VARIABLES
# 2 L1: OBJECT-ORIENTED IMAGE ANALYSIS - PART 2
#   ... Neighbor Operations: Class 44, 55, 66, and growing of class 111 (scarps)
#   ... Statistics of final output: Convexity, Flow, Main Direction 
#   ... Final Cleaning and Rasterisation: creation of class 11




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 1 PACKAGES, FUNCTIONS & VARIABLES ------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, maptools, cleangeo, mclust, raster, rgdal, rgeos, rgrass7, RSAGA, tools, sf, sp, shapefiles)

# path
setwd("E:/Masterarbeit/Data")




print("... Load Data")
# load("L1_part1.RData")
load("L1_part2.RData")



# ... initialize SAGA & GRASS --------------------------------------------------------
print("Initialize SAGA and GRASS")

# SAGA
env <-   rsaga.env(path = "../Software/saga_2.2.2_x64")
env$version # "2.2.2"
# env

# GRASS
epsg.code<-'EPSG:31256'
grass.epsg.code <- as.numeric(substring(epsg.code,6,10))  # define projection | MGI_Austria_GK_East
grass.loc <- paste0('loc',substring(epsg.code,6,10))      # define corresponding folder name
grass.mapset <- 'THESIS'   
grass.gis.base <- "/usr/local/grass-7.2.0" # server


# if grass location is not exist create a new grass location with the needed projection
if (!file.exists(file.path(getwd(), "Output/GrassData",grass.loc)))
{
  system(paste0("grass72 -e -c EPSG:", grass.epsg.code, " ", getwd(), "/Output/GrassData/" , grass.loc))
} 


# if mapset is not exisiting create and initialize new mapset
if (!file.exists(file.path(getwd(), "Output/GrassData",grass.loc, grass.mapset)))
{
  initGRASS(gisBase = grass.gis.base, home = tempdir(), 
            gisDbase = file.path(getwd(), "Output/GrassData"), override=TRUE, 
            location = grass.loc, mapset = grass.mapset, SG = as(raster::init(dtm.tif), "SpatialGrid"))
  
  
  # print(parseGRASS("r.in.gdal"))
  execGRASS('r.in.gdal',  flags=c('o',"overwrite", "quiet"), input=dtm.GRASS,  output='dtm', Sys_show.output.on.console = FALSE)
  
  
  # print(parseGRASS("g.region"))
  execGRASS("g.region",  raster="dtm", Sys_show.output.on.console = FALSE)
} else
{
  initGRASS(gisBase = grass.gis.base, home = tempdir(), 
            gisDbase = file.path(getwd(), "Output/GrassData"), override=TRUE, 
            location = grass.loc, mapset = grass.mapset, SG = as(raster::init(dtm.tif), "SpatialGrid"))
}

gmeta()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ... FUNCTIONS ------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
source("../R/server_module_Segmentation.R") # segmentation.tmp.path is temp folder directory
source("../R/server_module_ObjectOrientation.R") 
source("../R/module_NeighborOperations.R")

kMeanThresholds <- function(df, n.clust = NULL, G = 25, size.data = NULL, size.sample = 5000, seed = 123, iter.max = 5000)
{
  
  if(any(is.na(df)) == TRUE)
  {
    if(class(df) == "data.frame" || class(df)[1] == "data.table")
    {
      df <- df[complete.cases(df),]
    } else
    {
      df <- df[complete.cases(df)]
    }
    
  }
  
  m <- as.matrix(df)
  
  
  if(!is.null(size.data))
  {
    if(!findInterval(size.data, c(1, 99) ) == 1)
    {
      stop("Wrong input for size.data. Only numbers in the range of 1 - 99 are allowed!")
    }
    
    set.seed(seed)
    s <- sample(1:nrow(m), size = (nrow(m) * size.data/100))
    m <-  as.matrix(m[s, ])
  }
  
  
  if(is.null(n.clust))
  {
    if(!is.null(size.sample) && nrow(m) >= size.sample)
    {
      set.seed(seed)
      #        User      System verstrichen 
      #       277.97        1.23      281.60 
      n.clust <- mclust::Mclust(m, G = 1:G, initialization = list(subset = sample(1:nrow(m), size = size.sample)))
    } 
    
    if(is.null(size.sample) || nrow(m) < size.sample)
    {
      n.clust <- mclust::Mclust(m, G = 1:G,  na.action=na.exclude)
    }
    
    n.clust.best <- dim(n.clust$z)[2]
  } else
  {
    n.clust.best <- n.clust
  }
  
  kmeans.cluster  <- kmeans(x = m, centers = n.clust.best, iter.max = iter.max)
  
  if( kmeans.cluster$ifault == 4)
  {
    kmeans.cluster <- kmeans(x = m, centers = n.clust.best, algorithm = "MacQueen", iter.max = iter.max) 
  }
  
  
  return(kmeans.cluster$centers)
}






correctDBF <- function(x, end.n = length(colnames(d$dbf)), adjust.n = 0, new.colnames)
{
  d <- suppressMessages(shapefiles::read.dbf(paste0(file_path_sans_ext(x), ".dbf")))
  colnames(d$dbf) <- c(colnames(d$dbf)[0:(end.n-adjust.n)], new.colnames)
  shapefiles::write.dbf(d, paste0(file_path_sans_ext(x), ".dbf")) # write dbf with better header
  rm(d)
}



replaceInvalids <- function(x, replace.value = -9999) 
{
  # # NA and NAN
  # x[is.na(x)] <- replace.value
  # 
  # # NULL
  # x[is.null(x)] <- replace.value
  
  # found here: http://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table
  # NA and NAN
  for(j in seq_len(ncol(x)))
    set(x, which(is.na(x[[j]])), j, replace.value)

  # NULL
  for(j in seq_len(ncol(x)))
    set(x, which(is.null(x[[j]])), j, replace.value)
  
  return(x)
}











# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 2 OBJECT-ORIENTED IMAGE ANALYSIS - PART 2 ------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# get thresholds

l1.thrs.dtm51dif[order(l1.thrs.dtm51dif),]
# -5.4716001 -2.5721536 -0.8806598  0.1519328  0.8253774  2.2043472

l1.thrs.entr[order(l1.thrs.entr),]
l1.thrs.sa[order(l1.thrs.sa),]

l1.thrs.slp[order(l1.thrs.slp),] # 23.314519
# 3.977664  8.010430 12.316246 17.242991 23.314519 30.767841 39.882225 

l1.thrs.srough[order(l1.thrs.srough),]
# 0.001026526 0.004200564 0.008403560 0.013420473 0.019348083 0.026205967 
# 0.033922443 0.042618731 0.053989406 0.068193191 0.085388775 0.107442800 0.137431652 0.184856559  


l1.thrs.svf[order(l1.thrs.svf),]
# 0.6446079 0.7432842 0.8167058 0.8744564 0.9196185 0.9578650 

l1.thrs.curv15 <- c(-0.0123849036, -0.0065062248, -0.0030500411, -0.0009440737,  0.0003991937,  0.0014039410,  0.0023107739,  0.0031985304,
                    0.0041145528,  0.0051595209,  0.0063694280,  0.0077523111,  0.0093509959,  0.0111553422,  0.0132056692,  0.0155216383,
                    0.0181607517,  0.0211461151,  0.0246543204,  0.0288220734,  0.0341879664,  0.0417965644,  0.0585871504)

l1.thrs.curv15[order(l1.thrs.curv15),]
# -0.0123849036 -0.0065062248 -0.0030500411 -0.0009440737  0.0003991937  0.0014039410  0.0023107739  0.0031985304 
# 0.0041145528  0.0051595209  0.0063694280  0.0077523111  0.0093509959  0.0111553422  0.0132056692  0.0155216383 
# 0.0181607517  0.0211461151  0.0246543204  0.0288220734  0.0341879664  0.0417965644  0.0585871504 

rm(seg1.sel.sf)
seg1.sel.sf <- sf::st_read(dsn = paste(getwd(), dirname(seg1.sel), sep = "/"), 
                           layer = file_path_sans_ext(basename(seg1.sel)), quiet = TRUE)





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# L1: Refine Data - Class 111 ----------------------------------------
print("... First Selection of possible Scarp Candidates - Class 111")
# "Area" >= 25 AND  "Slp" >= 23.314519 AND "DTM_D_51" >= -2.5721536 AND "DTM_D_51" <= 0.8253774 AND "SVF" >= 0.7432842 AND  "SVF" <= 0.8744564
seg1.sel.sf$Class[seg1.sel.sf$Area >= 25 & seg1.sel.sf$Slp  >=  23.314519 & seg1.sel.sf$DTM_D_51  >= -2.5721536  & seg1.sel.sf$DTM_D_51  < 0.8253774 & seg1.sel.sf$SVF >= 0.7432842  & seg1.sel.sf$SVF < 0.8744564 & seg1.sel.sf$SfRghn >= 0.004200564] <- 111


# subset by class
# seg1.sel.sf.sub <- subset(seg1.sel.sf, seg1.sel.sf$Class == 111)

# convert sf to sp for further analysis
seg1.sel.sp <- as(seg1.sel.sf, "Spatial")







# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# L1: Neighbor Operation ----------------------------------------

# # #  Get Geometric Information

# geometry of Class 111
L1.class111 <- subset(seg1.sel.sp, seg1.sel.sp@data$Class == 111)


# ... SVF and Oppenness: Class 44 ----------------------------------------

# get intersecting IDs from SVF and Oppenness
# subsetting and union data
L1.SVF.sub <- subset(seg1.sel.sp, seg1.sel.sp@data$SVF < 0.7432842)
# L1.SVF.sub.union <- rgeos::gUnaryUnion(L1.SVF.sub)
# L1.SVF.sub.union  <- checkGeometry(L1.SVF.sub.union)

L1.dtmD51.sub <- subset(seg1.sel.sp, seg1.sel.sp@data$DTM_D_51 < -2.5721536)
# L1.dtmD51.sub.union <- rgeos::gUnaryUnion(L1.dtmD51.sub)
# L1.dtmD51.sub.union  <- checkGeometry(L1.dtmD51.sub.union)

# intersect geometries
# L1.SVF.inters.dtmD51 <- unlist(gIntersects(spgeom1 = L1.SVF.sub.union, spgeom2 = L1.dtmD51.sub, byid = TRUE, returnDense = FALSE))
# L1.dtmD51.inters.SVF <- unlist(gIntersects(spgeom1 = L1.dtmD51.sub.union, spgeom2 = L1.SVF.sub, byid = TRUE, returnDense = FALSE))
L1.SVF.inters.dtmD51 <- unlist(gIntersects(spgeom1 = L1.SVF.sub, spgeom2 = L1.dtmD51.sub, byid = TRUE, returnDense = FALSE))
L1.dtmD51.inters.SVF <- unlist(gIntersects(spgeom1 = L1.dtmD51.sub, spgeom2 = L1.SVF.sub, byid = TRUE, returnDense = FALSE))


# get ID's out of dtmD51 subset
L1.class44.ID <- unique(c(L1.dtmD51.sub[L1.SVF.inters.dtmD51,]@data$ID, L1.SVF.sub[L1.dtmD51.inters.SVF,]@data$ID))
L1.class44.ID.pos <- match(L1.class44.ID, seg1.sel.sp@data$ID)


# intersection to class44
L1.class44 <- seg1.sel.sp[L1.class44.ID.pos,]
# L1.class44.union <- rgeos::gUnaryUnion(L1.class44)
# L1.class44.union  <- checkGeometry(L1.class44.union)

# L1.class44.inters.class111 <- unlist(gIntersects(spgeom1 = L1.class44.union, spgeom2 = L1.class111, byid = TRUE, returnDense = FALSE))
L1.class44.inters.class111 <- unlist(gIntersects(spgeom1 = L1.class44, spgeom2 = L1.class111, byid = TRUE, returnDense = FALSE))


# get ID's out of L1.class111 subset and get position afterwards
L1.class44.inters.class111.ID <- L1.class111[L1.class44.inters.class111,]@data$ID
L1.class44.inters.class111.pos <- match(L1.class44.inters.class111.ID, seg1.sel.sp@data$ID)

# set class 44 into class 111
seg1.sel.sp@data$Class[L1.class44.inters.class111.pos] <- 44



# ... Curvature: Class 55 ----------------------------------------

# get intersecting IDs from SVF and Oppenness
# subsetting and union data
L1.curv.sub <- subset(seg1.sel.sp, seg1.sel.sp@data$CurvM15 >=  0.00411455)

# intersection to class55 (without union faster!)
L1.class55.inters.class111 <- unlist(gIntersects(spgeom1 = L1.curv.sub, spgeom2 = L1.class111, byid = TRUE, returnDense = FALSE))


# get ID's out of L1.class111 subset and get position afterwards
L1.class55.inters.class111.ID <- L1.class111[-L1.class55.inters.class111,]@data$ID # - because the interest lay on objects without border to high maximal curvature values
L1.class55.inters.class111.pos <- match(L1.class55.inters.class111.ID, seg1.sel.sp@data$ID)

# set class 55 into class 111
seg1.sel.sp@data$Class[L1.class55.inters.class111.pos] <- 55



# rgdal::writeOGR(subset(seg1.sel.sp, seg1.sel.sp@data$Class == 111 | seg1.sel.sp@data$Class == 44 | seg1.sel.sp@data$Class == 55 ), dsn = paste(getwd(), dirname(seg1.sel), sep = "/"), overwrite_layer = TRUE,
#                 layer = "L1_class111_4455", driver = "ESRI Shapefile")



# ... Stream: Class 66 ----------------------------------------
# intersection to class66 (without union faster!)
# read class66
stream <- "Output/Other/flow_Dinf_stream.shp"
L1.class66 <- sf::st_read(dsn = paste(getwd(), dirname(stream), sep = "/"), 
                           layer = file_path_sans_ext(basename(stream)), quiet = TRUE)

L1.class66 <- as(L1.class66, "Spatial")

L1.class66.inters.class111 <- unlist(gIntersects(spgeom1 = L1.class66, spgeom2 = L1.class111, byid = TRUE, returnDense = FALSE))


# get ID's out of L1.class111 subset and get position afterwards
L1.class66.inters.class111.ID <- L1.class111[L1.class66.inters.class111,]@data$ID
L1.class66.inters.class111.pos <- match(L1.class66.inters.class111.ID, seg1.sel.sp@data$ID)


# set class 66 into class 111
seg1.sel.sp@data$Class[L1.class66.inters.class111.pos] <- 66


# rgdal::writeOGR(subset(seg1.sel.sp, seg1.sel.sp@data$Class == 111 | seg1.sel.sp@data$Class == 44 | seg1.sel.sp@data$Class == 55 | seg1.sel.sp@data$Class == 66), dsn = paste(getwd(), dirname(seg1.sel), sep = "/"), overwrite_layer = TRUE,
#                 layer = "L1_class111_445566", driver = "ESRI Shapefile")





# ... Neighbor Growing of Class 111: Candidates of Scarp Parts to single polygon ----------------------------------------

# create neihgborhood
# L1_nb_speed_up <- rgeos::gUnarySTRtreeQuery(seg1.sel.sp) # speed up function for poly2nb
# L1_nb <- spdep::poly2nb(seg1.sel.sp, queen = TRUE, foundInBox = nb_speed_up) # neighborhood based on queen continuity


# grow scarp parts togehter and write them out
# 20 Min
L1.NeighborGrowing <- NeighborGrowing(spdf = seg1.sel.sp, ID.class = subset(seg1.sel.sp@data, seg1.sel.sp@data$Class == 111)$ID,
                                      return.input = TRUE, return.gUnaryUnionNeighbors = TRUE)


# necessairy to write growing data!
rgdal::writeOGR(L1.NeighborGrowing[[2]], dsn = paste(getwd(), dirname(L1.final), sep = "/"), overwrite_layer = TRUE, 
                layer = file_path_sans_ext(basename(L1.final)), driver = "ESRI Shapefile")




# # # # # # # # # 
# seg1.sel.sp.out <- seg1.sel.sp
# replaceInvalids(seg1.sel.sp.out@data)
# 
# rgdal::writeOGR(seg1.sel.sp.out, dsn = paste(getwd(), dirname(seg1.sel), sep = "/"), overwrite_layer = TRUE,
#                 layer = "L1_seg_sel_out", driver = "ESRI Shapefile")
# # # # # # # # # # 




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# L1: Statistics of final output: Convexity, Flow, Main Direction ----------------------------------------


L1.final.CH <- "Output/Segmentation/L1_final_ch.shp"
L1.final.pts <- "Output/Segmentation/L1_final_pts.shp"

# shape indices
rsaga.geoprocessor(lib="shapes_polygons", module = 7, env = env, show.output.on.console = FALSE, param = list(
  SHAPES = L1.final, INDEX = L1.final))


# compute convex hull and their shape indices
# convex hull
# rsaga.get.usage("shapes_points", 12, env = env)
rsaga.geoprocessor(lib="shapes_points", module = 12, env = env, show.output.on.console = FALSE, param = list(
  SHAPES = L1.final, HULLS = L1.final.CH))

rsaga.geoprocessor(lib="shapes_polygons", module = 7, env = env, show.output.on.console = FALSE, param = list(
  SHAPES = L1.final.CH, INDEX = L1.final.CH))


# flow direction
# flow.Dinf.cos <- paste0(file_path_sans_ext(flow.Dinf.cos), ".sgrd")
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = flow.Dinf.cos, POLYGONS = L1.final, METHOD = "0", NAMING = "1", RESULT = L1.final,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))

# flow.Dinf.sin <- paste0(file_path_sans_ext(flow.Dinf.sin), ".sgrd")
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = flow.Dinf.sin, POLYGONS = L1.final, METHOD = "0", NAMING = "1", RESULT = L1.final,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))





# ... shape indices and convex hull ----------------------------------------

# read shapes and convert them into sp
L1.final.CH.sf <- sf::st_read(dsn = paste(getwd(), dirname(L1.final.CH), sep = "/"), layer = file_path_sans_ext(basename(L1.final.CH)), quiet = TRUE, stringsAsFactors = FALSE)
L1.final.CH.sp <- as(L1.final.CH.sf, 'Spatial')

L1.final.sf <- sf::st_read(dsn = paste(getwd(), dirname(L1.final), sep = "/"), layer = file_path_sans_ext(basename(L1.final)), quiet = TRUE, stringsAsFactors = FALSE)
L1.final.sp <- as(L1.final.sf, 'Spatial')



# calculate convexity and compactness
L1.final.sp@data$Conv <- L1.final.CH.sp@data$Area/L1.final.CH.sp@data$Area.1
L1.final.sp@data$Comp <- L1.final.sp@data$Area/(L1.final.sp@data$Perimeter^2) *4 * pi



colnames(L1.final.sp@data) <- c("ID", "Area", "P","P_A", "P_sqrt_A", "Mx_Dist", "D_A", "D_sqrt_A", "Sh_Ind", 
                                "Fl_Cos", "Fl_Sin", "Conv", "Comp")






#  ... main direction and flow direction ----------------------------------------
print("... Computation of main flow and mean direction")

# # # mean flow
# set NAs
L1.final.sp@data$Fl_Cos[((L1.final.sp@data$Fl_Cos == 0) & (L1.final.sp@data$Fl_Sin == 0)) | ((L1.final.sp@data$Fl_Cos == -9999) & (L1.final.sp@data$Fl_Sin == -9999))] <- NA
L1.final.sp@data$Fl_Sin[is.na(L1.final.sp@data$Fl_Cos)] <- NA

# calculate mean Flow
L1.final.sp@data$Flow <- ((atan2(L1.final.sp@data$Fl_Sin, L1.final.sp@data$Fl_Cos) * (-180)/pi) + 90) %% 360



# # # mean flow: caluclation of object orientation and objects lying perpendicular to flow direction
L1.final.MainDir <- MainDirection(L1.final.sp)
L1.final.sp@data$MnDir <- L1.final.MainDir$angle
L1.final.sp@data$MnDirInv <- L1.final.MainDir$angle.inv

L1.final.sp@data$MnFl_D <- abs(L1.final.sp@data$MnDir - L1.final.sp@data$Flow)
L1.final.sp@data$MnInfFl_D <- abs(L1.final.sp@data$MnDirInv - L1.final.sp@data$Flow)

percentage.FlowMinPer <- 30
L1.final.sp@data$FlMnPer<- ifelse(((L1.final.sp@data$MnFl_D >= (90*((100-percentage.FlowMinPer)/100)) & L1.final.sp@data$MnFl_D <= (90*((100+percentage.FlowMinPer)/100)) ) | (L1.final.sp@data$MnInfFl_D >= (90*((100-percentage.FlowMinPer)/100)) & L1.final.sp@data$MnInfFl_D <= (90*((100+percentage.FlowMinPer)/100)))), 1, 0)


# calculation of length-width-ratio
L1.final.LenWidthRatio <- LengthWidthRatio(L1.final.sp)
L1.final.sp@data$LeWiRat <- L1.final.LenWidthRatio$ratio


# assign flow perpendicluar to class 11
L1.final.sp@data$Class <- NA
# L1.final.sp@data$Class[which(L1.final.sp@data$FlMnPer == 1)] <- 11







# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# L1: Final Cleaning and Rasterisation ------------------------------------
# final cleaning


#  ... get thresholds of shape characteristics ----------------------------------------
# convexity
l1.thrs.conv <- kMeanThresholds(L1.final.sp@data$Conv, size.sample = NULL)
l1.thrs.conv[order(l1.thrs.conv),]
# 0.4667959 0.6868825 
# 0.3081035 0.4939739 0.6280104 0.7493972 


# compactness
l1.thrs.comp <- kMeanThresholds(L1.final.sp@data$Comp, size.sample = NULL)
l1.thrs.comp[order(l1.thrs.comp),]
# 0.07409421 0.14756019 0.23155630 0.34751494 
# 0.07212148 0.14708015 0.23074610 0.35023983 


# shape index
l1.thrs.Sh_Ind <- kMeanThresholds(L1.final.sp@data$Sh_Ind, size.sample = NULL)
l1.thrs.Sh_Ind[order(l1.thrs.Sh_Ind),]
# 2.082387 2.902504 4.057335 6.112338 
# 2.002056 2.669966 3.491007 4.692490 6.930348 


# D_A
l1.thrs.D_A <- kMeanThresholds(L1.final.sp@data$D_A, size.sample = NULL)
l1.thrs.D_A[order(l1.thrs.D_A),]
# 0.1274540 0.2941322 0.4696749 
# 0.1261285 0.2859841 0.4588673 


# LenWiRat
l1.thrs.LeWiRat <- kMeanThresholds(L1.final.sp@data$LeWiRat, size.sample = NULL)
l1.thrs.LeWiRat[order(l1.thrs.LeWiRat),]
# 6.612209   23.232190   47.635474   86.942898  145.223219  220.831783  336.999609  581.615315 1319.960052 
#  6.298345   20.049473   39.420032   66.574313  101.485930  155.262427  220.657163  335.747764  566.713725 1319.960052  

# Area
#l1.thrs.Area <- kMeanThresholds(L1.final.sp@data$Area, size.sample = NULL)
#l1.thrs.Area[order(l1.thrs.Area),]
# 84.38471   301.01282   615.47632  1090.40437  1783.44330  2834.79032  4457.39024  7515.86364 13360.20000 
# 76.67937   252.38970   505.40378   897.52888  1416.88136  2111.40206  3297.48864  5232.21212  7587.47826 13131.38462 





#  ... selection of class 11 as final of L1 ----------------------------------------

# OLD: "FlMnPer" =1 AND "LeWiRat" < 95.063132  AND "Area" < 5241.85294  AND NOT ("Comp" <= 0.0733 AND "Conv" > 0.4898145) 
# L1.final.sp@data$Class[which(L1.final.sp@data$FlMnPer == 1 & L1.final.sp@data$LeWiRat < 145.223219 & L1.final.sp@data$Area < 4457.39024 & !(L1.final.sp@data$Comp < 0.07409421 & L1.final.sp@data$Conv >= 0.4667959))] <- 11

# "FlMnPer" = 1 AND "LeWiRat" < 145.223219  AND "Area" < 5000 AND "Sh_Ind" < 4.692490 OR ("FlMnPer" = 0 AND "Conv" < 0.4667959  AND "Comp"  >= 0.07409421) 
L1.final.sp@data$Class[which(L1.final.sp@data$FlMnPer == 1 & L1.final.sp@data$LeWiRat < 145.223219 & L1.final.sp@data$Sh_Ind < 4.7 & L1.final.sp@data$Area < 5000 | (L1.final.sp@data$FlMnPer == 0 & L1.final.sp@data$Conv < 0.4667959 & L1.final.sp@data$Comp >= 0.07409421))] <- 11



# # # # # # # # # 
# L1.final.sp.out <- L1.final.sp
# replaceInvalids(L1.final.sp.out@data)
# 
# rgdal::writeOGR(L1.final.sp.out, dsn = paste(getwd(), dirname(seg1.sel), sep = "/"), overwrite_layer = TRUE,
#                 layer = "L1_final_out", driver = "ESRI Shapefile")
# # # # # # # # # # 


source("../R/module_NeighborOperations.R")

# subset to class 11
L1.final.sp.sub <- subset(L1.final.sp, L1.final.sp@data$Class == 11)

# correct projection
L1.final.sp.sub <- spTransform(L1.final.sp.sub, CRS(paste0("+init=", epsg.code))) 


# # correct geometry
L1.final.sp.sub <- checkGeometry(L1.final.sp.sub)


# write out data
replaceInvalids(L1.final.sp.sub@data, replace.value = -9999)
rgdal::writeOGR(L1.final.sp.sub, dsn = paste(getwd(), dirname(L1.final), sep = "/"), overwrite_layer = TRUE, 
                layer = file_path_sans_ext(basename(L1.final)), driver = "ESRI Shapefile")




#  ... creation of points in L1_final and buffering L1_final for next segmentation step ----------------------------------------

# re-read shapefile to get clear geometry
L1.final.mapT <- readShapePoly(fn = L1.final, proj4string = CRS(paste0("+init=", epsg.code)))

# create points and its dataframe
L1.final.points.sp <- rgeos::gPointOnSurface(L1.final.mapT, byid = TRUE)
L1.final.points.spdf <- SpatialPointsDataFrame(L1.final.points.sp, data.frame(index = row.names(L1.final.mapT)))

# write out points
rgdal::writeOGR(L1.final.points.spdf, dsn = paste(getwd(), dirname(L1.final.pts), sep = "/"), overwrite_layer = TRUE, 
                layer = file_path_sans_ext(basename(L1.final.pts)), driver = "ESRI Shapefile")


# rasterize shapefile
# rsaga.get.usage("grid_gridding", 0, env = env)
# POLY_TYPE:  [0] node
# OUTPUT: [1] index number
# GRID_TYPE: [2] Integer (4 byte)
# TARGET_DEFINITION: grid or grid system
rsaga.geoprocessor(lib="grid_gridding", module = 0, env = env, show.output.on.console = FALSE, param = list(
  INPUT = L1.final, OUTPUT = "1", POLY_TYPE = "0", GRID_TYPE = "2", TARGET_DEFINITION = "1", TARGET_TEMPLATE = dtm, GRID = L1.final.grid))



if(!file.exists(file.path(getwd(), L1.final.grid.buf)))
{
  # grid buffering
  # rsaga.get.usage("grid_tools", 8, env = env)
  # BUFFERTYPE: [0] Fixed | [1] Cell value 
  rsaga.geoprocessor(lib="grid_tools", module = 8, env = env, show.output.on.console = FALSE, param = list(
    FEATURES = L1.final.grid, BUFFER = L1.final.grid.buf, DIST = L1.final.grid.buf.dist, BUFFERTYPE = "0")) 

}


# change storage format
# rsaga.get.usage("grid_tools", 11, env = env)
# TYPE:  [5] unsigned 4 byte integer
# rsaga.geoprocessor(lib="grid_tools", module = 11, env = env, show.output.on.console = FALSE, param = list(
#   INPUT = L1.final.grid, OUTPUT = L1.final.grid, TYPE = "5"))


# redefine no data to zero
# rsaga.get.usage("grid_calculus", 1, env = env)
# TYPE:[5] unsigned 4 byte integer
rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = FALSE, param=list(
  GRIDS = L1.final.grid, RESULT = L1.final.grid, FORMULA = "ifelse(g1= (-99999), 0, g1)",
  FNAME = "1", USE_NODATA = "1", TYPE = "5"))





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Save Data ---------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
print("... Save Data")
save.image("L1_part2.RData")






