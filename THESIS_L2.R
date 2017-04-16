# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# THESIS: SECOND SEGMENTATION LEVEL ------------------------------------------------------------------------
#
# Raphael Knevels
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# DESCRIPTION:
# In this script the object-oriented image analysis for segmentation level 2
# is performed. The objective of this levels is the classification of potential
# flank objects. Moreover, crowns are extracted.
#
# This script builds on the final outputs of in THESIS_L1_part2.R.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #




# CONTENT -----------------------------------------------------------------
# 1 PACKAGES, FUNCTIONS & VARIABLES
# 2 L2: OBJECT-ORIENTED IMAGE ANALYSIS
#   ... fine-scale segmentation level 2
#   ... pre-selection based on L1 results
#   ... declaring of variables of L2
#   ... calculate grid statistics openness, maximal curvature and shape indices 
#   ... Neighbor and Class Operations: flank-scarp and crown-scarp relation
#   ... calculation of k-means thresholds
#   ... Selection of Flanks and Crowns and Modifying: class 22 and 77
#   ... Neighbor Growing, Points and Rasterisation: Create final output in L2




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
load("L2.RData")



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
source("../R/module_NeighborOperations.R")
source("../R/module_ObjectOrientation.R")
source("../R/module_BoundingBox.R")


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

}











# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# OBJECT-ORIENTED IMAGE ANALYSIS ------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ... SEGMENTATION LEVEL 2 ------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

print("Start SEGMENTATION LEVEL 2 ...")

Segments.Grid.2 <- "Output/Segmentation/Sec2/seg2_01_15_real.sgrd"
Segments.Poly.2 <- "Output/Segmentation/Sec2/seg2_01_15_real.shp"
Segments.Poly.2.sel <- "Output/Segmentation/Sec2/seg2_sel_01_15_real.shp"
segmentation(Tool = "GRASS", Segments.Grid = Segments.Grid.2, Segments.Poly = Segments.Poly.2, Input.Grid = c(dtm.dif.dtm51, slope.edge.vigra, "<>", dtm.dif.dtm51, slope.edge.vigra.re03),
             Seed.Method = "Fast Representativeness", Fast.Representativeness.LevelOfGeneralisation = 1.5,  Saga.Segmentation.Method = "0", Saga.Segmentation.Sig.2 = "125",
             Grass.Segmentation.Minsize = 25, Grass.Segmentation.Threshold = "0.001", Grass.Segmentation.Memory = 3072, Segmentation.Boundary.Grid = L1.final.grid, Grass.Segmentation.Weighted = TRUE,
             Generalisation.Flac = FALSE, Saga.Segmentation.Leafsize = 1024, NoData = TRUE, Mask = slope, burn.Boundary.into.Segments = c(TRUE, TRUE))





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ... pre-selection based on L1 results ------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = L1.final.grid.buf, POLYGONS = Segments.Poly.2, METHOD = "0", NAMING = "1", RESULT = Segments.Poly.2,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "1", MEAN = "0", VAR = "0", STDDEV = "0", QUANTILE = 0))





# correct corrupt field naming of SAGA GIS
correctDBF(Segments.Poly.2, adjust.n = 1, new.colnames = c("L1_Sum"))


seg2.sf <- sf::st_read(dsn = paste(getwd(), dirname(Segments.Poly.2), sep = "/"), layer = file_path_sans_ext(basename(Segments.Poly.2)), quiet = TRUE, stringsAsFactors = FALSE)
seg2.sel.sp <- as(seg2.sel.sf, "Spatial")

# select all polygon that intersected with the max curvature raster
seg2.sel.sp <- seg2.sel.sp[seg2.sel.sp@data$L1_Sum > 0,]

# write to shapefile
rgdal::writeOGR(seg2.sel.sp, dsn = paste(getwd(), dirname(Segments.Poly.2.sel), sep = "/"), overwrite_layer = TRUE, 
                layer = file_path_sans_ext(basename(Segments.Poly.2.sel)), driver = "ESRI Shapefile")







# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ... declaring of variables of L2  ------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# 
slope <- "Output/Slope/slope.sgrd"
dtm <- "Output/DTM/dtm.sgrd"
dtm.dif.dtm51 <- "Output/DTM/dtm_dif_dtm51.sgrd"
curv.max15 <- "Output/Other/curv_max15.sgrd"
flow.Dinf.cos <- "Output/Other/flow_Dinf_cos.tif"
flow.Dinf.sin <-  "Output/Other/flow_Dinf_sin.tif"
epsg.code<-'EPSG:31256'

# # # new variables
L1.final <- "Output/Segmentation/L1_final.shp"
L1.final.pts <- "Output/Segmentation/L1_final_pts.shp"
L1.final.pts.sf  <- sf::st_read(dsn = paste(getwd(), dirname(L1.final.pts), sep = "/"), layer = file_path_sans_ext(basename(L1.final.pts)), quiet = TRUE, stringsAsFactors = FALSE)
L1.final.pts.sp <- as(L1.final.pts.sf, "Spatial") 
L1.final.grid <- "Output/Segmentation/L1_final.sgrd"


L1.final.sf <- sf::st_read(dsn = paste(getwd(), dirname(L1.final), sep = "/"), layer = file_path_sans_ext(basename(L1.final)), quiet = TRUE, stringsAsFactors = FALSE)
L1.final.sp <- as(L1.final.sf, "Spatial") 
L1.final.tbl <- L1.final.sp@data



L2.flank <- "Output/Segmentation/L2_flank.shp"
L2.flank.pts <- "Output/Segmentation/L2_flank_pts.shp"
L2.flank.grid <- "Output/Segmentation/L2_flank.sgrd"

L2.crown <- "Output/Segmentation/L2_crown.shp"
L2.crown.pts <- "Output/Segmentation/L2_crown_pts.shp"
L2.crown.grid <- "Output/Segmentation/L2_crown.sgrd"

L2.final.grid <- "Output/Segmentation/L2_final.sgrd"
L2.final.grid.buf <- "Output/Segmentation/L2_final_buf.sgrd"
L2.final.grid.buf.dist <- 25


# convert sf to sp for further analysis
seg2.sel <- "Output/Segmentation/L2_seg_sel.shp"



# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# L2: calculate grid statistics dtm dif 51, maximal curvature and shape indices ----------------------------------------
# rsaga.get.usage("shapes_grid", 2, env = env)


# slope
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = slope, POLYGONS = seg2.sel, METHOD = "0", NAMING = "1", RESULT = seg2.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))



# dtm dif 51
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = dtm.dif.dtm51, POLYGONS = seg2.sel, METHOD = "0", NAMING = "1", RESULT = seg2.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))



# maximal curvature 
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = curv.max15, POLYGONS = seg2.sel, METHOD = "0", NAMING = "1", RESULT = seg2.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# flow direction
# flow.Dinf.cos <- paste0(file_path_sans_ext(flow.Dinf.cos), ".sgrd")
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = flow.Dinf.cos, POLYGONS = seg2.sel, METHOD = "0", NAMING = "1", RESULT = seg2.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))

# flow.Dinf.sin <- paste0(file_path_sans_ext(flow.Dinf.sin), ".sgrd")
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = flow.Dinf.sin, POLYGONS = seg2.sel, METHOD = "0", NAMING = "1", RESULT = seg2.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))



# shape indices
rsaga.geoprocessor(lib="shapes_polygons", module = 7, env = env, show.output.on.console = FALSE, param = list(
  SHAPES = seg2.sel, INDEX = seg2.sel))



correctDBF(seg2.sel, adjust.n = 13, new.colnames = c("Slp","DTM_D_51", "CurvM15", "Fl_Cos", "Fl_Sin", "Area", "P","P_A", "P_sqrt_A", "Mx_Dist", "D_A", "D_sqrt_A", "Sh_Ind"))




# ... main direction, length-width-ratio and flow
# caluclation of object orientation and objects lying perpendicular to flow direction

seg2.sel.sf <- sf::st_read(dsn = paste(getwd(), dirname(seg2.sel), sep = "/"), layer = file_path_sans_ext(basename(seg2.sel)), quiet = TRUE, stringsAsFactors = FALSE)
seg2.sel.sp <- as(seg2.sel.sf, "Spatial")

# main direction
# "------ Run of MainDirection: 2.03250000000001 Minutes ------"
seg2.sel.MainDir <- MainDirection(seg2.sel.sp)
seg2.sel.sp@data$MnDir <- seg2.sel.MainDir$angle
seg2.sel.sp@data$MnDirInv <- seg2.sel.MainDir$angle.inv


# length width ratio
seg2.sel.sp@data$LeWiRat <- LengthWidthRatio(seg2.sel.sp)$ratio


# caluclation of Flow Direction and Inverse
seg2.sel.sp@data$Fl_Cos[((seg2.sel.sp@data$Fl_Cos == 0) & (seg2.sel.sp@data$Fl_Sin == 0)) | ((seg2.sel.sp@data$Fl_Cos == -9999) & (seg2.sel.sp@data$Fl_Sin == -9999))] <- NA
seg2.sel.sp@data$Fl_Sin[is.na(seg2.sel.sp@data$Fl_Cos)] <- NA
# seg1.sel.sf$Fl_Sin[seg1.sel.sf$Fl_Cos == -9999] <- -9999
seg2.sel.sp@data$Flow <- ((atan2(seg2.sel.sp@data$Fl_Sin, seg2.sel.sp@data$Fl_Cos) * (-180)/pi) + 90) %% 360
seg2.sel.sp@data$FlowInv <-(seg2.sel.sp@data$Flow  - 180 + 360) %% 360 



# calculate difference between flow and main direction
seg2.sel.sp@data$MnFl_D <- abs(seg2.sel.sp@data$MnDir - seg2.sel.sp@data$Flow)
seg2.sel.sp@data$MnInfFl_D <- abs(seg2.sel.sp@data$MnDirInv - seg2.sel.sp@data$Flow)


# get class 1 from L1 in L2 (main scarps)
L2.class1.pos <- unlist(rgeos::gIntersects(L1.final.pts.sp,  seg2.sel.sp, byid = TRUE, returnDense = FALSE))
seg2.sel.sp@data$Class[L2.class1.pos] <- 11
L1.final.tbl$L2_pos <- L2.class1.pos


# # #
# transform NA to -9999
# seg2.sel.sp.out <- seg2.sel.sp
# replaceInvalids(seg2.sel.sp.out@data, replace.value = -9999)
# 
# rgdal::writeOGR(seg2.sel.sp.out, dsn = paste(getwd(), dirname(seg2.sel), sep = "/"), overwrite_layer = TRUE, 
#                 layer = "L2_seg_sel_out", driver = "ESRI Shapefile")
# # #


# L2: Neighbor and Class Operations --------------------------------------------

L2.nb.speed.up <- rgeos::gUnarySTRtreeQuery(seg2.sel.sp) # speed up function for poly2nb
L2.nb <- spdep::poly2nb(seg2.sel.sp, queen = TRUE, foundInBox = L2.nb.speed.up) # neighborhood based on queen continuity



# ... get flank neighbor relation ---------------------------------
# get relation of neighbors to class 11: main direction analysis to detect flank that are perpendicular to it (of interest: mean absolute difference to class)
# [1] "------ Run of RelationalClassFunction: 18.2071666666667 Minutes ------"
L2.MnDir.RelClass <- RelationalClassFunction(spdf = seg2.sel.sp, nb =  L2_nb, class.var = 11, var = c("MnDir", "Area"), var.is.angle = TRUE, class.as.neighbors = FALSE)



# ... get neighbors in flow and inverse flow direction ---------------------------------

# geometry of Class 11
# L2.class11 <- subset(seg2.sel.sp, seg2.sel.sp@data$Class == 11)


# calculation of bounding box in flow and inverse flow direction
L2.bb.flow <- getBoundingBox(shape = subset(seg2.sel.sp, seg2.sel.sp@data$Class == 11), col.name = "Flow", scale.factor = c(1.4, 0.4), k.centroid = 10,
                              centroid = FALSE, set.centroid = "direction", scale.side = "small", quiet = TRUE)


L2.bb.flowInv <- getBoundingBox(shape = subset(seg2.sel.sp, seg2.sel.sp@data$Class == 11), col.name = "FlowInv", scale.factor = 0.6, k.centroid = 4,
                              centroid = FALSE, set.centroid = "direction", scale.side = "small", quiet = TRUE)

# writeSpatialShape(x = L2.bb.flow, fn = paste(getwd(), dirname(seg2.sel), "/Second Segmentation/bb/bb_flow8", sep = "/"))
# writeSpatialShape(x = L2_bb_flowInv, fn = paste(getwd(), dirname(seg2.sel), "/Second Segmentation/bb/bb_flowInv", sep = "/"))


# get all objects that intersect the bounding boxes
L2.L2bbflow.inters.seg2 <- unique(unlist(gIntersects(spgeom1 = L2.bb.flow, spgeom2 = seg2.sel.sp, byid = TRUE, returnDense = FALSE)))
L2.L2bbflowInv.inters.seg2 <- unique(unlist(gIntersects(spgeom1 = L2.bb.flowInv, spgeom2 = seg2.sel.sp, byid = TRUE, returnDense = FALSE)))


# # kick out class 11 objects and objects that are not neighbors of class 11
L2.class11.pos <- which(seg2.sel.sp@data$Class == 11)
L2.nb.class11 <- unique(unlist(L2.nb[L2.class11.pos]))

# kick out class 11 objects
L2.L2bbflow.inters.seg2 <- L2.L2bbflow.inters.seg2[!L2.L2bbflow.inters.seg2 %in% L2.class11.pos]
L2.L2bbflowInv.inters.seg2 <- L2.L2bbflowInv.inters.seg2[!L2.L2bbflowInv.inters.seg2 %in% L2.class11.pos]

# kick out not neighbor objects
L2.L2bbflow.inters.seg2  <- L2.L2bbflow.inters.seg2[L2.L2bbflow.inters.seg2 %in% L2.nb.class11]
L2.L2bbflowInv.inters.seg2  <- L2.L2bbflowInv.inters.seg2[L2.L2bbflowInv.inters.seg2 %in% L2.nb.class11]

# put information to data 
seg2.sel.sp$NbCl11Fl <- 0
seg2.sel.sp$NbCl11FlInv <- 0

seg2.sel.sp$NbCl11Fl[L2.L2bbflow.inters.seg2] <- 1
seg2.sel.sp$NbCl11FlInv[L2.L2bbflowInv.inters.seg2] <- 1



# merge tables to SPDF
seg2.sel.sp@data <- merge(seg2.sel.sp@data, L2.MnDir.RelClass[c("ID", "ang_d_cl_a", "ang_m_cl")], by = "ID", suffixes = c("", "_MnD"), all.x = TRUE)
seg2.sel.sp@data$angInv_d_cl_a <- abs(seg2.sel.sp@data$MnDirInv - seg2.sel.sp@data$ang_m_cl)




# ... calculate flanks and crowns by thresholds ----------------------------------------
# ... calculation for flanks ---------------------------------
percentage.FlankMinPer <- 30
percentage.FlankMinPerConv <- 50
differenceFlankMainFlow <- 27

flank.conv.L1.pos <- L1.final.tbl[L1.final.tbl$Conv < 0.4667959,]$L2_pos
flank.conv.pos  <- unique(unlist(L2.nb[flank.conv.L1.pos]))

# Flank to Main Direction of class should be 90°
seg2.sel.sp@data$Flank <- ifelse(((seg2.sel.sp@data$ang_d_cl_a >= (90*((100-percentage.FlankMinPer)/100)) & seg2.sel.sp@data$ang_d_cl_a  <= (90*((100+percentage.FlankMinPer)/100)) ) | (seg2.sel.sp@data$angInv_d_cl_a >= (90*((100-percentage.FlankMinPer)/100)) & seg2.sel.sp@data$angInv_d_cl_a <= (90*((100+percentage.FlankMinPer)/100)))), 1, 0)
seg2.sel.sp@data$Flank[flank.conv.pos] <- ifelse(((seg2.sel.sp@data$ang_d_cl_a[flank.conv.pos] >= (90*((100-percentage.FlankMinPerConv)/100)) & seg2.sel.sp@data$ang_d_cl_a[flank.conv.pos]  <= (90*((100+percentage.FlankMinPerConv)/100)) ) | (seg2.sel.sp@data$angInv_d_cl_a[flank.conv.pos] >= (90*((100-percentage.FlankMinPerConv)/100)) & seg2.sel.sp@data$angInv_d_cl_a[flank.conv.pos] <= (90*((100+percentage.FlankMinPerConv)/100)))), 1, 0)


# Flank should not have same direction as flow
seg2.sel.sp@data$Flank[which(seg2.sel.sp@data$Flank == 1)] <- ifelse(((seg2.sel.sp@data$MnFl_D[which(seg2.sel.sp@data$Flank == 1)] >= differenceFlankMainFlow) | (seg2.sel.sp@data$MnInfFl_D[which(seg2.sel.sp@data$Flank == 1)] >= (180+differenceFlankMainFlow))), 1, 0)



# ... crowns ------------------------------
differenceCrownMain <- 27
seg2.sel.sp@data$Crown<- ifelse(((seg2.sel.sp@data$ang_d_cl_a  <= differenceCrownMain) | (seg2.sel.sp@data$ang_d_cl_a  >= (180-differenceCrownMain))  | (seg2.sel.sp@data$angInv_d_cl_a >= (360-differenceCrownMain)) | (seg2.sel.sp@data$angInv_d_cl_a <= (180+differenceCrownMain))), 1, 0)



# # # 
# seg2.sel.sp.out <- seg2.sel.sp
# replaceInvalids(seg2.sel.sp.out@data, replace.value = -9999)
# 
# rgdal::writeOGR(seg2.sel.sp.out, dsn = paste(getwd(), dirname(seg2.sel), sep = "/"), overwrite_layer = TRUE,
#                 layer = "L2_seg_sel_out", driver = "ESRI Shapefile")
# 
# correctDBF("Output/Segmentation/L2_seg_sel_out.shp", adjust.n = 28, new.colnames = c("Slp", "DTM_D_51", "CurvM15", "Fl_Cos", "Fl_Sin",
#                                                      "Area", "P","P_A", "P_sqrt_A", "Mx_Dist",
#                                                      "D_A", "D_sqrt_A", "Sh_Ind","MnDir", "MnDirInv", "LeWiRat",
#                                                      "Flow", "FlowInv", "MnFl_D", "MnInfFl_D" , "Class", "NbCl11Fl", "NbCl11FlInv",
#                                                      "angle", "angle_mCl","angle_inv", "Flank", "Crown"))
# # #



# # #
# L2: calculate k-means thresholds ----------------------------------------
# # #

# dtm to dtm51 dif
l2.thrs.dtm51dif <- kMeanThresholds(seg2.sel.sf$DTM_D_51)
l2.thrs.dtm51dif[order(l2.thrs.dtm51dif),]
# -1.8480002 -0.1382933  1.0989332 



# slp
# l1.thrs.dtm51dif <- kMeanThresholds(seg1.sel.sf[seg1.sel.sf$t_Entr >= l1.thrs.entr.sel,]$DTM51_Dif)
l2.thrs.slp <- kMeanThresholds(seg2.sel.sf$Slp)
l2.thrs.slp[order(l2.thrs.slp),]
# 7.072568 13.407496 19.588355 26.283062 34.615099 



# D_A
l2.thrs.D_A <- kMeanThresholds(seg2.sel.sf$D_A)
l2.thrs.D_A[order(l2.thrs.D_A),]
# 0.06748803 0.10131988 0.13805622 0.18113174 0.22809245 0.27978090 0.33752436 0.40357238 0.48833247 0.59517102 0.74288104 1.23358660 



# Shp_Index
l2.thrs.Sh_Ind <- kMeanThresholds(seg2.sel.sf$Sh_Ind)
l2.thrs.Sh_Ind[order(l2.thrs.Sh_Ind),]
# 1.925149 2.621215 3.453631 4.858446 


# CurvM15
options(scipen=999)
l2.thrs.CurvM15 <- kMeanThresholds(seg2.sel.sf$CurvM15[seg2.sel.sf$CurvM15 >= -1])
l2.thrs.CurvM15[order(l2.thrs.CurvM15),]
# 0.0006365585 0.0053492877 0.0115429194 0.0214774322 

# LeWiRat
l2.thrs.LeWiRat <- kMeanThresholds(seg2.sel.sp@data$LeWiRat)
l2.thrs.LeWiRat[order(l2.thrs.LeWiRat),]
#  3.027028    7.448578   13.321072   20.964592   30.692215   42.759480   58.601485   78.964692  103.913564  135.413909 176.960206  232.954370  315.449842  420.016491  578.089453  865.326689 1246.279782 







# # # 
# L2: Selection of Flanks and Crowns and Modifying ------------------------------------
# # #

# ... class 22: Flanks --------------------------------
# "Flank" = 1 AND "Slp" >= 13.407496 AND "DTM_D_51" >= -1.8480002 AND "Area" >= 25 AND "D_A" >=  0.10131988  AND "Sh_Ind" <  3.453631 AND "LeWiRat" < 58.601485 AND "CurvM15" < 0.0115429194 AND "NbCl11Fl" = 1
seg2.sel.sp@data$Class[seg2.sel.sp@data$Flank == 1 & seg2.sel.sp@data$NbCl11Fl == 1 & seg2.sel.sp@data$Area >= 25 & seg2.sel.sp@data$Slp  >=  13.407496 & seg2.sel.sp@data$DTM_D_51  >= -1.8480002  & seg2.sel.sp@data$D_A >= 0.10131988 & seg2.sel.sp@data$Sh_Ind < 3.453631 & seg2.sel.sp@data$LeWiRat < 58.601485 & seg2.sel.sp@data$CurvM15 < 0.0115429194 & seg2.sel.sp@data$Class!= 11] <- 22


# remove flanks that have an higher area then the scarp
# get relation of class 22 to 11: only flanks that have a smaller mean area are retained
# [1] "------ Run of RelationalClassFunction: 1.14350000000001 Minutes ------"
class22.to.11.area <- RelationalClassFunction(spdf = seg2.sel.sp, nb =  L2.nb, class.var = 11, class.var2 = 22,
                                            var = c("Area", "Area"), class.as.neighbors = FALSE, quiet = FALSE)

# merge tables to SPDF
seg2.sel.sp@data <- merge(seg2.sel.sp@data, class22.to.11.area[c("ID", "m_d_nb_cl")], by = "ID", suffixes = c("", "_A"), all.x = TRUE)


# kick out class 22 elements that has a higher area then the scarp (class 11)
seg2.sel.sp@data$Class[seg2.sel.sp@data$Class == 22 & seg2.sel.sp@data$m_d_nb_cl > 0] <- -9999




# ... class 77: Crowns --------------------------------------
# "Crown" = 1 AND "NbCl11FlIn" = 1 AND "Slp" < 19.588355 AND "CurvM15"  >= 0.0115429194
# seg2.sel.sp@data$Class[seg2.sel.sp@data$Crown == 1 & seg2.sel.sp@data$NbCl11FlInv == 1 & seg2.sel.sp@data$Area >= 25 & seg2.sel.sp@data$Slp  < 19.588355 &  seg2.sel.sp@data$CurvM15 >= 0.0115429194 & seg2.sel.sp@data$Class!= 11 & seg2.sel.sp@data$Class!= 22] <- 77
seg2.sel.sp@data$Class[seg2.sel.sp@data$NbCl11FlInv == 1 & seg2.sel.sp@data$Area >= 25 & seg2.sel.sp@data$Slp  < 19.588355 &  seg2.sel.sp@data$CurvM15 >= 0.0115429194 & seg2.sel.sp@data$Class!= 11 & seg2.sel.sp@data$Class!= 22] <- 77



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# L2: Neighbor Growing, Points and Rasterisation ------------------------------------

# ... Neighbor Growing --------------------------------------------
# growing of flanks (class 22)
L2.NeighborGrowing.cl22 <- NeighborGrowing(spdf = seg2.sel.sp, ID.class = subset(seg2.sel.sp@data, seg2.sel.sp@data$Class == 22)$ID, return.gUnaryUnionNeighbors = TRUE, return.input = FALSE)
                                      

# necessairy to write growing data!
rgdal::writeOGR(L2.NeighborGrowing.cl22, dsn = paste(getwd(), dirname(L2.flank), sep = "/"), overwrite_layer = TRUE, 
                layer = file_path_sans_ext(basename(L2.flank)), driver = "ESRI Shapefile")





# growing of crowns (class 77)
L2.NeighborGrowing.cl77 <- NeighborGrowing(spdf = seg2.sel.sp, ID.class = subset(seg2.sel.sp@data, seg2.sel.sp@data$Class == 77)$ID, return.gUnaryUnionNeighbors = TRUE, return.input = FALSE)


# necessairy to write growing data!
rgdal::writeOGR(L2.NeighborGrowing.cl77, dsn = paste(getwd(), dirname(L2.crown), sep = "/"), overwrite_layer = TRUE, 
                layer = file_path_sans_ext(basename(L2.crown)), driver = "ESRI Shapefile")





# ... get Points of Classes for the next segmentation steps  --------------------------------------------

# re-read shapefile to get clear geometry
L2.flank.mapT <- readShapePoly(fn = L2.flank, proj4string = CRS(paste0("+init=", epsg.code)))
L2.crown.mapT <- readShapePoly(fn = L2.crown, proj4string = CRS(paste0("+init=", epsg.code)))


# create points and its dataframe
L2.flank.pts.sp <- rgeos::gPointOnSurface(L2.flank.mapT, byid = TRUE)
L2.flank.pts.spdf <- SpatialPointsDataFrame(L2.flank.pts.sp, data.frame(index = row.names(L2.flank.mapT)))

L2.crown.pts.sp <- rgeos::gPointOnSurface(L2.crown.mapT, byid = TRUE)
L2.crown.pts.spdf <- SpatialPointsDataFrame(L2.crown.pts.sp, data.frame(index = row.names(L2.crown.mapT)))


# write out points
rgdal::writeOGR(L2.flank.pts.spdf, dsn = paste(getwd(), dirname(L2.flank.pts), sep = "/"), overwrite_layer = TRUE, 
                layer = file_path_sans_ext(basename(L2.flank.pts)), driver = "ESRI Shapefile")

rgdal::writeOGR(L2.crown.pts.spdf, dsn = paste(getwd(), dirname(L2.crown.pts), sep = "/"), overwrite_layer = TRUE, 
                layer = file_path_sans_ext(basename(L2.crown.pts)), driver = "ESRI Shapefile")



# ... Rasterisation  --------------------------------------------


# ... rasterize shapefile -------------------
# rsaga.get.usage("grid_gridding", 0, env = env)
# POLY_TYPE:  [0] node
# OUTPUT: [1] index number
# GRID_TYPE: [2] Integer (4 byte)
# TARGET_DEFINITION: grid or grid system
rsaga.geoprocessor(lib="grid_gridding", module = 0, env = env, show.output.on.console = FALSE, param = list(
  INPUT = L2.flank, OUTPUT = "1", POLY_TYPE = "0", GRID_TYPE = "2", TARGET_DEFINITION = "1", TARGET_TEMPLATE = dtm, GRID = L2.flank.grid))


rsaga.geoprocessor(lib="grid_gridding", module = 0, env = env, show.output.on.console = FALSE, param = list(
  INPUT = L2.crown, OUTPUT = "1", POLY_TYPE = "0", GRID_TYPE = "2", TARGET_DEFINITION = "1", TARGET_TEMPLATE = dtm, GRID = L2.crown.grid))




# ... redefine no data to zero and add values to keep unique IDs -------------------------
# rsaga.get.usage("grid_calculus", 1, env = env)
# TYPE:[5] unsigned 4 byte integer

L2.flank.val  <- nrow(L1.final.sf) + 8
L2.flank.formula.expression <- paste0("ifelse(g1 = 0, 0,  (g1+", L2.flank.val, "))")


rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = FALSE, param=list(
  GRIDS = L2.flank.grid, RESULT = L2.flank.grid, FORMULA = L2.flank.formula.expression,
  FNAME = "1", USE_NODATA = "1", TYPE = "5"))



L2.crown.val  <- nrow(L1.final.sf) + 8 + length(which(seg2.sel.sp@data$Class == 22)) + 8
L2.crown.formula.expression <- paste0("ifelse(g1 = 0, 0,  (g1+", L2.crown.val, "))")

rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = FALSE, param=list(
  GRIDS = L2.crown.grid, RESULT = L2.crown.grid, FORMULA = L2.crown.formula.expression,
  FNAME = "1", USE_NODATA = "1", TYPE = "5"))




# ... create L2 final output by adding all grids together

rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = FALSE, param=list(
  GRIDS = paste(c(L1.final.grid, L2.flank.grid, L2.crown.grid), collapse = ";"), RESULT = L2.final.grid, FORMULA = "a + b + c",
  FNAME = "1", USE_NODATA = "1", TYPE = "5"))




if(!file.exists(file.path(getwd(), L2.final.grid.buf)))
{
  # grid buffering
  # rsaga.get.usage("grid_tools", 8, env = env)
  # BUFFERTYPE: [0] Fixed | [1] Cell value 
  rsaga.geoprocessor(lib="grid_tools", module = 8, env = env, show.output.on.console = FALSE, param = list(
    FEATURES = L2.final.grid, BUFFER = L2.final.grid.buf, DIST = L2.final.grid.buf.dist, BUFFERTYPE = "0")) 
  
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Save Data ---------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
print("... Save Data")
save.image("L2.RData")






