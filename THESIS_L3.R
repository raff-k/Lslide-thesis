# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# THESIS: Third SEGMENTATION LEVEL ------------------------------------------------------------------------
#
# Raphael Knevels
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# DESCRIPTION:
# In this script the object-oriented image analysis for segmentation level 3
# is performed. The objective of this levels is the classification of potential
# body objects. Moreover, all landslide parts are further refinded, and finally
# grown into one single landslide.
#
# This script builds on the final outputs of in THESIS_L1_part2.R and THESIS_L2.R.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #




# CONTENT -----------------------------------------------------------------
# 1 PACKAGES, FUNCTIONS & VARIABLES
# 2 L3: OBJECT-ORIENTED IMAGE ANALYSIS
#   ... declaration of variables
#   ... coarser scale segmentation L3
#   ... pre-selection of L3 segmentation
#   ... calculate grid statistics: land surface parameters, textural features and shape indices
#   ... get classes out of previous segmentation levels: class 11, 22, 77
#   ... Neighbor and Class Operations: neighbors in flow direction, class 88, 89 and 99
#   ... calculation of k-means thresholds
#   ... refine class 11 to class 10
#   ... Landslide Body for Class 10: class 33 and class 1 - final scarps
#   ... Landslide Body - Class 3 and class 1, class 2 and class 7
#   ... Neighbor Growing and final Cleaning: Growing of class 1, 2 and 3 to final single landslide
# 3 Accuracy Assessment: graphs








# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 1 PACKAGES, FUNCTIONS & VARIABLES ------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, maptools, cleangeo, mclust, raster, rgdal, rgeos, rgrass7, ggplot2, RSAGA, tools, sf, sp, shapefiles)

# path
setwd("E:/Masterarbeit/Data")




print("... Load Data")
# load("L1_part1.RData")
load("L3.RData")



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
grass.gis.base <-  "../Software/GRASS GIS 7.2.0" # server


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
source("../R/module_ObjectFeatures.R")

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
# 2 OBJECT-ORIENTED IMAGE ANALYSIS ------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ... variables ------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

slope <- "Output/Slope/slope.sgrd"
slope.texture <- "Output/Slope/slope_texture.tif"
dtm <- "Output/DTM/dtm.sgrd"
dtm.dif.dtm51 <- "Output/DTM/dtm_dif_dtm51.sgrd"
skyviewfactor <- "Output/DTM/dtm_SVF_R20_D16.tif"
surf.roughn15 <- "Output/Other/surf_roughn15.sgrd"
curv.max15 <- "Output/Other/curv_max15.sgrd"

flow.Dinf.deg <- "Output/Other/flow_Dinf_deg.tif" 
texture.var  <- "Output/Other/textureFlowDir51_VAR.tiff" 
texture.contr <- "Output/Other/textureFlowDir51_Contr.tiff" 
texture.entr <- "Output/Other/textureFlowDir51_Entr.tiff" 
texture.idm <- "Output/Other/textureFlowDir51_IDM.tiff" 

TextureFlow(grass.texture.flowDir = paste0(getwd(), "/", flow.Dinf.deg), grass.texture.input = "slope", grass.texture.method = "contrast,var,idm,entr", grass.texture.name = c("_Contr_", "_VAR_", "_IDM_","_Entr_"),
            grass.texture.window = 5, grass.texture.distance = 1, grass.texture.save = TRUE, grass.texture.save.path = texture.entr.flow)


# # # normalized height
dtm.r3 <- "Output/DTM/dtm_r3.sgrd"
dtm.r3.NH <- "Output/DTM/dtm_r3_NH_5_2_2.sgrd"

# rsaga.get.usage("grid_tools", 0, env = env)
# SCALE_UP: [3] Bicubic Spline Interpolation
rsaga.geoprocessor(lib="grid_tools", module = 0, env = env, show.output.on.console = FALSE, param = list(
  INPUT = dtm, OUTPUT = dtm.r3, SCALE_UP = "3", TARGET_DEFINITION = "0", TARGET_USER_SIZE = "3"))

# rsaga.get.usage("ta_morphometry", 14, env = env)
rsaga.geoprocessor(lib="ta_morphometry", module = 14, env = env, show.output.on.console = FALSE, param = list(
  DEM = dtm.r3, NH = dtm.r3.NH, W = "5", T = "2", E = "2"))




flow.Dinf.cos <- "Output/Other/flow_Dinf_cos.tif"
flow.Dinf.sin <-  "Output/Other/flow_Dinf_sin.tif"
epsg.code<-'EPSG:31256'

# # # L1 and L2 variables
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
L2.flank.pts.sf  <- sf::st_read(dsn = paste(getwd(), dirname(L2.flank.pts), sep = "/"), layer = file_path_sans_ext(basename(L2.flank.pts)), quiet = TRUE, stringsAsFactors = FALSE)
L2.flank.pts.sp <- as(L2.flank.pts.sf, "Spatial") 
L2.flank.grid <- "Output/Segmentation/L2_flank.sgrd"

L2.crown <- "Output/Segmentation/L2_crown.shp"
L2.crown.pts <- "Output/Segmentation/L2_crown_pts.shp"
L2.crown.pts.sf  <- sf::st_read(dsn = paste(getwd(), dirname(L2.crown.pts), sep = "/"), layer = file_path_sans_ext(basename(L2.crown.pts)), quiet = TRUE, stringsAsFactors = FALSE)
L2.crown.pts.sp <- as(L2.crown.pts.sf, "Spatial") 
L2.crown.grid <- "Output/Segmentation/L2_crown.sgrd"

L2.final.grid <- "Output/Segmentation/L2_final.sgrd"

# convert sf to sp for further analysis
seg2.sel <- "Output/Segmentation/L2_seg_sel.shp"


# # # new L3 variables
seg3.grid <- "Output/Segmentation/L3_seg.sgrd"
seg3.sel <- "Output/Segmentation/L3_seg_sel.shp"
seg3.sel.CH <- "Output/Segmentation/L3_seg_sel_CH.shp"

seg3.sel.grid <- "Output/Segmentation/L3_seg_sel.tif"

L3.landslide <- "Output/Segmentation/L3_landslide.shp"





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ... SEGMENTATION LEVEL 2 ------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


Segments.Grid.2 <- "Output/Segmentation/L3_seg.sgrd"
Segments.Poly.2 <- "Output/Segmentation/L3_seg.shp"
Segments.Poly.2.sel <- "Output/Segmentation/L3_seg_sel.shp"

segmentation(Tool = "GRASS", Segments.Grid = Segments.Grid.3, Segments.Poly = Segments.Poly.3, Input.Grid = c(slope, slope.edge.vigra, "<>", slope, slope.edge.vigra.re03),
            Seed.Method = "Fast Representativeness", Fast.Representativeness.LevelOfGeneralisation = 4.25,  Saga.Segmentation.Method = "0", Saga.Segmentation.Sig.2 = "125",
            Grass.Segmentation.Minsize = 25, Grass.Segmentation.Threshold = "0.0001", Grass.Segmentation.Memory = 10000, Segmentation.Boundary.Grid = L2.final.grid, Grass.Segmentation.Weighted = TRUE,
            Generalisation.Flac = FALSE, Saga.Segmentation.Leafsize = 1024, NoData = TRUE, Mask = slope, burn.Boundary.into.Segments = c(TRUE, TRUE))



rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = L2.final.grid.buf, POLYGONS = Segments.Poly.2, METHOD = "0", NAMING = "1", RESULT = Segments.Poly.2,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "1", MEAN = "0", VAR = "0", STDDEV = "0", QUANTILE = 0))





# correct corrupt field naming of SAGA GIS
correctDBF(Segments.Poly.3, adjust.n = 1, new.colnames = c("L2_Sum"))


# ... pre-selection of L3 segmentation ------------------------------------------
seg3.sf <- sf::st_read(dsn = paste(getwd(), dirname(Segments.Poly.3), sep = "/"), layer = file_path_sans_ext(basename(Segments.Poly.3)), quiet = TRUE, stringsAsFactors = FALSE)
seg3.sel.sp <- as(seg2.sf, "Spatial")

# select all polygon that intersected with the max curvature raster
seg3.sel.sp <- seg3.sel.sp[seg3.sel.sp@data$L2_Sum > 0,]

# write to shapefile
rgdal::writeOGR(seg3.sel.sp, dsn = paste(getwd(), dirname(Segments.Poly.3.sel), sep = "/"), overwrite_layer = TRUE, 
                layer = file_path_sans_ext(basename(Segments.Poly.3.sel)), driver = "ESRI Shapefile")









# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# L3: calculate grid statistics: land surface parameters, textural features and shape indices ----------------------------------------
# rsaga.get.usage("shapes_grid", 2, env = env)

# slope
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = slope, POLYGONS = seg3.sel, METHOD = "0", NAMING = "1", RESULT = seg3.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "1", QUANTILE = 0))


# CurvMax 
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = curv.max15, POLYGONS = seg3.sel, METHOD = "0", NAMING = "1", RESULT = seg3.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# dtm
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = dtm, POLYGONS = seg3.sel, METHOD = "0", NAMING = "1", RESULT = seg3.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# dtm dif 51
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = dtm.dif.dtm51, POLYGONS = seg3.sel, METHOD = "0", NAMING = "1", RESULT = seg3.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# normalized height
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = dtm.r3.NH, POLYGONS = seg3.sel, METHOD = "0", NAMING = "1", RESULT = seg3.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# flow direction
# flow.Dinf.cos <- paste0(file_path_sans_ext(flow.Dinf.cos), ".sgrd")
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = flow.Dinf.cos, POLYGONS = seg3.sel, METHOD = "0", NAMING = "1", RESULT = seg3.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))

# flow.Dinf.sin <- paste0(file_path_sans_ext(flow.Dinf.sin), ".sgrd")
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = flow.Dinf.sin, POLYGONS = seg3.sel, METHOD = "0", NAMING = "1", RESULT = seg3.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))



# vector ruggedness measure 
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = surf.roughn15, POLYGONS = seg3.sel, METHOD = "0", NAMING = "1", RESULT = seg3.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))




# sky view factor
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = skyviewfactor, POLYGONS = seg3.sel, METHOD = "0", NAMING = "1", RESULT = seg3.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))



# texture entropy
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = texture.entr, POLYGONS = seg3.sel, METHOD = "0", NAMING = "1", RESULT = seg3.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# texture idm
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = texture.idm, POLYGONS = seg3.sel, METHOD = "0", NAMING = "1", RESULT = seg3.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# texture contrast
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = texture.contr, POLYGONS = seg3.sel, METHOD = "0", NAMING = "1", RESULT = seg3.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# texture variance
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = texture.var, POLYGONS = seg3.sel, METHOD = "0", NAMING = "1", RESULT = seg3.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# shape indices
rsaga.geoprocessor(lib="shapes_polygons", module = 7, env = env, show.output.on.console = FALSE, param = list(
  SHAPES = seg3.sel, INDEX = seg3.sel))


# convex hull and their shape indices
# convex hull
# rsaga.get.usage("shapes_points", 12, env = env)
rsaga.geoprocessor(lib="shapes_points", module = 12, env = env, show.output.on.console = FALSE, param = list(
  SHAPES = seg3.sel, HULLS = seg3.sel.CH))

rsaga.geoprocessor(lib="shapes_polygons", module = 7, env = env, show.output.on.console = FALSE, param = list(
  SHAPES = seg3.sel.CH, INDEX = seg3.sel.CH))



correctDBF(seg3.sel, adjust.n = 22, new.colnames = c("Slp", "Slp_SD", "CurvM15", "DTM", "DTM_D_51",  "DTM_NH", "Fl_Cos", "Fl_Sin", "SfRghn", "SVF",
                                                     "fl_Entr", "fl_IDM","fl_Contr", "fl_VAR", "Area", "P","P_A", "P_sqrt_A", "Mx_Dist", 
                                                     "D_A", "D_sqrt_A", "Sh_Ind"))






# ... main direction, length-width-ratio and flow ----------------------
# caluclation of object orientation and objects lying perpendicular to flow direction

seg3.sel.sf <- sf::st_read(dsn = paste(getwd(), dirname(seg3.sel), sep = "/"), layer = file_path_sans_ext(basename(seg3.sel)), quiet = TRUE, stringsAsFactors = FALSE)
seg3.sel.sp <- as(seg3.sel.sf, "Spatial")

# main direction
# [1] "------ Run of MainDirection: 2.6235 Minutes ------"
seg3.sel.MainDir <- MainDirection(seg3.sel.sp)
seg3.sel.sp@data$MnDir <- seg3.sel.MainDir$angle
seg3.sel.sp@data$MnDirInv <- seg3.sel.MainDir$angle.inv


# length width ratio
# [1] "------ Run of LengthWidthRatio: 1.56016666666668 Minutes ------"
seg3.sel.LeWiRat <- LengthWidthRatio(seg3.sel.sp)
seg3.sel.sp@data$LeWiRat <- seg3.sel.LeWiRat$ratio



# caluclation of Flow Direction and Inverse
seg3.sel.sp@data$Fl_Cos[((seg3.sel.sp@data$Fl_Cos == 0) & (seg3.sel.sp@data$Fl_Sin == 0)) | ((seg3.sel.sp@data$Fl_Cos == -9999) & (seg3.sel.sp@data$Fl_Sin == -9999))] <- NA
seg3.sel.sp@data$Fl_Sin[is.na(seg3.sel.sp@data$Fl_Cos)] <- NA

seg3.sel.sp@data$Flow <- ((atan2(seg3.sel.sp@data$Fl_Sin, seg3.sel.sp@data$Fl_Cos) * (-180)/pi) + 90) %% 360
seg3.sel.sp@data$FlowInv <-(seg3.sel.sp@data$Flow  - 180 + 360) %% 360 



# calculate difference between flow and main direction
seg3.sel.sp@data$MnFl_D <- abs(seg3.sel.sp@data$MnDir - seg3.sel.sp@data$Flow)
seg3.sel.sp@data$MnInfFl_D <- abs(seg3.sel.sp@data$MnDirInv - seg3.sel.sp@data$Flow)



# ... calculate compactness and convexity ---------------

# read shapes and convert them into sp
seg3.sel.CH.sf <- sf::st_read(dsn = paste(getwd(), dirname(seg3.sel.CH), sep = "/"), layer = file_path_sans_ext(basename(seg3.sel.CH)), quiet = TRUE, stringsAsFactors = FALSE)
seg3.sel.CH.sp <- as(seg3.sel.CH.sf, 'Spatial')


# calculate convexity and compactness
seg3.sel.sp@data$Conv <- seg3.sel.CH.sp@data$Area/seg3.sel.CH.sp@data$Area.1
seg3.sel.sp@data$Comp <- seg3.sel.sp@data$Area/(seg3.sel.sp@data$P^2) *4 * pi






# ... calculate haralick texture -------------------------
# # clip segmentation output raster first
# rsaga.get.usage("shapes_grid", 7, env = env)
rsaga.geoprocessor(lib="shapes_grid", module = 7, env = env, show.output.on.console = FALSE, param = list(
  OUTPUT = paste0(file_path_sans_ext(seg3.sel.grid), ".sgrd"), INPUT = seg3.grid, POLYGONS = seg3.sel))

# change format to tiff
# rsaga.get.usage("io_gdal", 1, env = env)
# FORMAT: [7] GeoTIFF
rsaga.geoprocessor(lib="io_gdal", module = 1, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = paste0(file_path_sans_ext(seg3.sel.grid), ".sgrd"), FILE = seg3.sel.grid, FORMAT = "7", SET_NODATA = "1", NODATA = "0"))

# ... same for slope
rsaga.geoprocessor(lib="io_gdal", module = 1, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = slope, FILE = slope.texture, FORMAT = "7", SET_NODATA = "1", NODATA = "0"))



# # calculation of haralick's texture
# [1] "------ Run of TextureObjectHaralick: 3.04533333333335 Minutes ------"
dt.haralick <- TextureObjectHaralick(texture.haralick.input = slope.texture, object.feature.segments = seg3.sel.grid, NA.val.in = 0, 
                                     texture.haralick.scales = c(1,5), texture.haralick.nbins = 32)


# merge to table
seg3.sel.sp@data <- merge(seg3.sel.sp@data, dt.haralick[, c("ID", "h_ent_s1", "h_idm_s1", "h_con_s5","h_var_s1"), with=FALSE], by.x = "L3_seg", by.y = "ID", all.x = TRUE)





# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# L3: get classes out of previous segmentation levels ----------------------------------------
# # # 

seg3.sel.sp@data$Class <- -9999


# get class 11 (main scarps)
L3.class11.pos <- unlist(rgeos::gIntersects(L1.final.pts.sp,  seg3.sel.sp, byid = TRUE, returnDense = FALSE))
seg3.sel.sp@data$Class[L3.class11.pos] <- 11


# get class 22 (flanks)
L3.class22.pos <- unlist(rgeos::gIntersects(L2.flank.pts.sp,  seg3.sel.sp, byid = TRUE, returnDense = FALSE))
seg3.sel.sp@data$Class[L3.class22.pos] <- 22


# get class 77 (crown)
L3.class77.pos <- unlist(rgeos::gIntersects(L2.crown.pts.sp,  seg3.sel.sp, byid = TRUE, returnDense = FALSE))
seg3.sel.sp@data$Class[L3.class77.pos] <- 77








# # #
# transform NA to -9999
# seg3.sel.sp.out <- seg3.sel.sp
# replaceInvalids(seg3.sel.sp.out@data, replace.value = -9999)
# 
# rgdal::writeOGR(seg3.sel.sp.out, dsn = paste(getwd(), dirname(seg3.sel), sep = "/"), overwrite_layer = TRUE,
#                 layer = "L3_seg_sel_out", driver = "ESRI Shapefile")
# # #






# # #
# L3: Neighbor and Class Operations --------------------------------------------
# # #


L3.nb.speed.up <- rgeos::gUnarySTRtreeQuery(seg3.sel.sp) # speed up function for poly2nb
L3.nb <- spdep::poly2nb(seg3.sel.sp, queen = TRUE, foundInBox = L3.nb.speed.up) # neighborhood based on queen continuity



# get class objects in L3
L3.class11.pos <- which(seg3.sel.sp@data$Class == 11)
L3.class22.pos <- which(seg3.sel.sp@data$Class == 22)
L3.class77.pos <- which(seg3.sel.sp@data$Class == 77)


# get neighbors of class objects in L3
L3.nb.class11 <- unique(unlist(L3.nb[L3.class11.pos]))
L3.nb.class22 <- unique(unlist(L3.nb[L3.class22.pos]))
L3.nb.class77 <- unique(unlist(L3.nb[L3.class77.pos]))




# ... get neighbors in flow direction ---------------------------------

# ... ... bounding box in flow direction ---------------------------------
# calculation of bounding box in flow and inverse flow direction
L3.bb.flow.Conv <- getBoundingBox(shape = subset(seg3.sel.sp, seg3.sel.sp@data$Conv < 0.4667959 & seg3.sel.sp@data$Class == 11), scale.factor = c(1.4, 1.2), 
                                  k.centroid = 3, set.centroid = "inverse", scale.side = "long", centroid = FALSE, col.name = "Flow", quiet = TRUE)

row.names(L3.bb.flow.Conv) <- row.names(subset(seg3.sel.sp, seg3.sel.sp@data$Conv < 0.4667959 & seg3.sel.sp@data$Class == 11))



L3.bb.flow.Centr <- getBoundingBox(shape = subset(seg3.sel.sp, seg3.sel.sp@data$Conv >= 0.4667959 & seg3.sel.sp@data$Class == 11), scale.factor = c(1.8, 1.2), 
                                   k = 2, scale.side = "long", centroid = TRUE, col.name = "Flow", quiet = TRUE)

row.names(L3.bb.flow.Centr) <- row.names(subset(seg3.sel.sp, seg3.sel.sp@data$Conv >= 0.4667959 & seg3.sel.sp@data$Class == 11))


# merge both bounding boxes
L3.bb.flow.cl11 <- rbind(L3.bb.flow.Conv, L3.bb.flow.Centr)

# order bounding box to fit ID classes order
L3.bb.flow.cl11  <- L3.bb.flow.cl11[order(L3.bb.flow.cl11$ID),]

# a <- subset(seg3.sel.sp, seg3.sel.sp@data$Class == 11)@data$ID
# b <- match(a, seg3.sel.sp@data$ID)

writeSpatialShape(x = L3.bb.flow.cl11, fn = paste(getwd(), dirname(seg3.sel), "L3_bb_flow_cl11.shp", sep = "/"))




# ... ... get neighbors in flow direction ---------------------------------
# get all objects that intersect the bounding boxes
L3.bbFlowCl11.inters.seg3 <- unique(unlist(gIntersects(spgeom1 = L3.bb.flow, spgeom2 = seg3.sel.sp, byid = TRUE, returnDense = FALSE)))


# kick out class objects
L3.bbFlowCl11.inters.seg3 <- L3.bbFlowCl11.inters.seg3[!L3.bbFlowCl11.inters.seg3 %in% c(L3.class11.pos, L3.class22.pos, L3.class77.pos)]


# kick out not neighbor objects
L3.bbFlowCl11.inters.seg3.nb  <- L3.bbFlowCl11.inters.seg3[L3.bbFlowCl11.inters.seg3 %in% c(L3.nb.class11)] # L3.nb.class22, L3.nb.class77



# put information to data 
seg3.sel.sp$Cl11Fl <- 0
seg3.sel.sp$Cl11FlNb <- 0

seg3.sel.sp$Cl11Fl[L3.bbFlowCl11.inters.seg3] <- 1
seg3.sel.sp$Cl11FlNb[L3.bbFlowCl11.inters.seg3.nb] <- 1






# ... neighbor operations ---------------------------------

# border of normalized height to class 11
# 0.08239538 0.19024493
# ... ... ... class 88 and 89: normalized height ---------------------------------
seg3.sel.sp@data$Class[which(seg3.sel.sp@data$Class == -9999 & seg3.sel.sp@data$DTM_NH < 0.08239538)] <- 88
L3.NH1.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb,  class.var = 88, class.var2 = 11, var = c("DTM_NH", "Area"), quiet = FALSE, class.as.neighbors = FALSE)

seg3.sel.sp@data$Class[which((seg3.sel.sp@data$Class == -9999 | seg3.sel.sp@data$Class == 88) & seg3.sel.sp@data$DTM_NH < 0.19024493)] <- 89
L3.NH2.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb,  class.var = 89, class.var2 = 11, var = c("DTM_NH", "Area"), quiet = FALSE, class.as.neighbors = TRUE)

L3.NH.RelBor.Cl11 <- cbind(L3.NH1.RelClass[c("ID", "Cl","bor_cl_rel")], L3.NH2.RelClass[c("bor_cl_rel")])
L3.NH.RelBor.Cl11 <- L3.NH.RelBor.Cl11[which(L3.NH.RelBor.Cl11$Cl == 11),]
L3.NH.RelBor.Cl11[is.na(L3.NH.RelBor.Cl11)] <- 0
colnames(L3.NH.RelBor.Cl11) <- c("ID", "Cl", "NH_rB_11_s", "NH_rB_11_l")
L3.NH.RelBor.Cl11$Cl <- NULL



# length(seg3.sel.sp@data$Class[which((seg3.sel.sp@data$Class == 88  | seg3.sel.sp@data$Class == 3 | seg3.sel.sp@data$Class == 89 | seg3.sel.sp@data$Class == -9999) & seg3.sel.sp@data$DTM_NH < 0.08239538)] )
# 2195

# length(seg3.sel.sp@data$Class[which((seg3.sel.sp@data$Class == 88  | seg3.sel.sp@data$Class == 3 | seg3.sel.sp@data$Class == 89 | seg3.sel.sp@data$Class == -9999) & seg3.sel.sp@data$DTM_NH < 0.19024493)] )
# 5836

# border of object with length to class 11
# 18.194977
# ... ... ... class 99: LeWiRat ---------------------------------
# [1] "------ Run of RelationalClassFunction: 2.07616666666666 Minutes ------"
seg3.sel.sp@data$Class[which(seg3.sel.sp@data$Class == -9999 & seg3.sel.sp@data$LeWiRat >= 38.942666)] <- 99
L3.LeWiRat.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb,  class.var = 99, class.var2 = 11, var = c("LeWiRat", "Area"), quiet = FALSE, class.as.neighbors = FALSE)

L3.LWR.RelBor.Cl11 <- L3.LeWiRat.RelClass[which(L3.LeWiRat.RelClass$Cl == 11),][c("ID", "bor_cl_rel")]
L3.LWR.RelBor.Cl11[is.na(L3.LWR.RelBor.Cl11)] <- 0
colnames(L3.LWR.RelBor.Cl11) <- c("ID", "LWR_rB_11")



# length(seg3.sel.sp@data$Class[which((seg3.sel.sp@data$Class == 99  | seg3.sel.sp@data$Class == 3 | seg3.sel.sp@data$Class == -9999) & seg3.sel.sp@data$LeWiRat >= 38.942666)] )
# 1102


# ... ... ClassNeighbor Function --------------------------------------------------
# statistics of neighbors to a given class


# surface roughness
# [1] "------ Run of ClassNeighborFunction: 1.10516666666666 Minutes ------"
L3.SfRghn.ClassNb <- ClassNeighborFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 11, var = c("SfRghn", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl11)
L3.SfRghn.ClassNb.11 <-  L3.SfRghn.ClassNb[c("ID", "m_wei", "m_w_FlMaN", "m_w_FlAlN", "m_dif", "m_dif_FlMa", "m_dif_FlAl")]
colnames(L3.SfRghn.ClassNb.11) <- c("ID", "SfRn_mw", "SfRn_mwfM", "SfRn_mwfA", "SfRn_md", "SfRn_mdfM", "SfRn_mdfA" )


# sky-view factor
L3.SVF.ClassNb <- ClassNeighborFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 11, var = c("SVF", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl11)
L3.SVF.ClassNb.11 <-  L3.SVF.ClassNb[c("ID", "m_wei", "m_w_FlMaN", "m_w_FlAlN", "m_dif", "m_dif_FlMa", "m_dif_FlAl")]
colnames(L3.SVF.ClassNb.11) <- c("ID", "SVF_mw", "SVF_mwfM", "SVF_mwfA", "SVF_md", "SVF_mdfM", "SVF_mdfA" )



# slope
L3.Slp.ClassNb <- ClassNeighborFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 11, var = c("Slp", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl11)
L3.Slp.ClassNb.11 <-  L3.Slp.ClassNb[c("ID", "m_wei", "m_w_FlMaN", "m_w_FlAlN", "m_dif", "m_dif_FlMa", "m_dif_FlAl", "sd_nb", "sd_nb_FlMaN" ,"sd_nb_FlAlN")]
colnames(L3.Slp.ClassNb.11) <- c("ID", "Slp_mw", "Slp_mwfM", "Slp_mwfA", "Slp_md", "Slp_mdfM", "Slp_mdfA", "Slp_sdcl", "Slp_sdfM", "Slp_sdfA")




# Openness
L3.dtmD51.ClassNb <- ClassNeighborFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 11, var = c("DTM_D_51", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl11)
L3.dtmD51.ClassNb.11 <-  L3.dtmD51.ClassNb[c("ID", "m_wei", "m_w_FlMaN", "m_w_FlAlN", "m_dif", "m_dif_FlMa", "m_dif_FlAl")]
colnames(L3.dtmD51.ClassNb.11) <- c("ID", "dtmD_mw", "dtmD_mwfM", "dtmD_mwfA", "dtmD_md", "dtmD_mdfM", "dtmD_mdfA")



# entropy flow
L3.flEntr.ClassNb <- ClassNeighborFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 11, var = c("fl_Entr", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl11)
L3.flEntr.ClassNb.11 <-  L3.flEntr.ClassNb[c("ID", "m_wei", "m_w_FlMaN", "m_w_FlAlN", "m_dif", "m_dif_FlMa", "m_dif_FlAl")]
colnames(L3.flEntr.ClassNb.11) <- c("ID", "fEnt_mw", "fEnt_mwfM", "fEnt_mwfA", "fEnt_md", "fEntr_mdfM", "fEnt_mdfA" )




# entropy object
L3.hEntr.ClassNb <- ClassNeighborFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 11, var = c("h_ent_s1", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl11)
L3.hEntr.ClassNb.11 <-  L3.hEntr.ClassNb[c("ID", "m_wei", "m_w_FlMaN", "m_w_FlAlN", "m_dif", "m_dif_FlMa", "m_dif_FlAl")]
colnames(L3.hEntr.ClassNb.11) <- c("ID", "hEnt_mw", "hEnt_mwfM", "hEnt_mwfA", "hEnt_md", "hEnt_mdfM", "hEnt_mdfA" )




# IDM flow
L3.flIDM.ClassNb <- ClassNeighborFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 11, var = c("fl_IDM", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl11)
L3.flIDM.ClassNb.11 <-  L3.flIDM.ClassNb[c("ID", "m_wei", "m_w_FlMaN", "m_w_FlAlN", "m_dif", "m_dif_FlMa", "m_dif_FlAl")]
colnames(L3.flIDM.ClassNb.11) <- c("ID", "fIDM_mw", "fIDM_mwfM", "fIDM_mwfA", "fIDM_md", "fIDM_mdfM", "fIDM_mdfA" )




# contrast flow
L3.flContr.ClassNb <- ClassNeighborFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 11, var = c("fl_Contr", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl11)
L3.flContr.ClassNb.11 <-  L3.flContr.ClassNb[c("ID", "m_wei", "m_w_FlMaN", "m_w_FlAlN", "m_dif", "m_dif_FlMa", "m_dif_FlAl")]
colnames(L3.flContr.ClassNb.11) <- c("ID", "fCont_mw", "fCont_mwfM", "fCont_mwfA", "fCont_md", "fCont_mdfM", "fCont_mdfA" )



# contrast object
L3.hContr.ClassNb <- ClassNeighborFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 11, var = c("h_con_s5", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl11)
L3.hContr.ClassNb.11 <-  L3.hContr.ClassNb[c("ID", "m_wei", "m_w_FlMaN", "m_w_FlAlN", "m_dif", "m_dif_FlMa", "m_dif_FlAl")]
colnames(L3.hContr.ClassNb.11) <- c("ID", "hCont_mw", "hCont_mwfM", "hCont_mwfA", "hCont_md", "hCont_mdfM", "hCont_mdfA" )



# variance flow
L3.flVar.ClassNb <- ClassNeighborFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 11, var = c("fl_VAR", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl11)
L3.flVar.ClassNb.11 <-  L3.flVar.ClassNb[c("ID", "m_wei", "m_w_FlMaN", "m_w_FlAlN", "m_dif", "m_dif_FlMa", "m_dif_FlAl")]
colnames(L3.flVar.ClassNb.11) <- c("ID", "fVar_mw", "fVar_mwfM", "fVar_mwfA", "fVar_md", "fVar_mdfM", "fVar_mdfA" )





# final merging fo tables, exclude ID columns
L3.ClassNb.11 <- cbind(L3.SfRghn.ClassNb.11, L3.SVF.ClassNb.11[-1], L3.Slp.ClassNb.11[-1], L3.dtmD51.ClassNb.11[-1], L3.flEntr.ClassNb.11[-1], 
                       L3.hEntr.ClassNb.11[-1], L3.flIDM.ClassNb.11[-1], L3.flContr.ClassNb.11[-1], L3.hContr.ClassNb.11[-1],
                       L3.flVar.ClassNb.11[-1])



# ...  merge to table -----------------------------------
seg3.sel.sp@data <- merge(seg3.sel.sp@data, L3.ClassNb.11, by = "ID", all.x = TRUE)
seg3.sel.sp@data <- merge(seg3.sel.sp@data, L3.NH.RelBor.Cl11, by = "ID", all.x = TRUE)
seg3.sel.sp@data <- merge(seg3.sel.sp@data, L3.LWR.RelBor.Cl11, by = "ID", all.x = TRUE)


# seg3.sel.sp@data <- seg3.sel.sp@data[-c( (ncol(seg3.sel.sp@data) - 62): ncol(seg3.sel.sp@data))]
# seg3.sel.sp@data <- seg3.sel.sp@data[-c( (ncol(seg3.sel.sp@data) - 2): ncol(seg3.sel.sp@data))]






# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# L3: calculate k-means thresholds ----------------------------------------
# # # 

# ... land surface parameter --------------------------------

# dtm51 dif
l3.thrs.dtm51dif <- kMeanThresholds(seg3.sel.sp@data$DTM_D_51)
l3.thrs.dtm51dif[order(l3.thrs.dtm51dif),]
# -3.2177916 -1.1489897  0.1280833  1.5434541
# -2.183391
# -3.2177916 -2.183391 -1.1489897  0.1280833  1.5434541


# normalized height
l3.thrs.dtmNH <- kMeanThresholds(seg3.sel.sp@data$DTM_NH[which(seg3.sel.sp@data$DTM_NH >= 0)])
l3.thrs.dtmNH[order(l3.thrs.dtmNH),]
# 0.08239538 0.19024493 0.29894804 0.41329221 0.53951430 0.68107097 0.84494128 



# slp
l3.thrs.slp <- kMeanThresholds(seg3.sel.sp@data$Slp)
l3.thrs.slp[order(l3.thrs.slp),]
#  4.606080  8.821035 13.321149 18.094487 23.176670 27.504224 32.481872 39.046244 



# Slp_SD 
l3.thrs.slp.sd <- kMeanThresholds(seg3.sel.sp@data$Slp_SD)
l3.thrs.slp.sd[order(l3.thrs.slp.sd),]
# 1.741864 3.311322 4.919841 6.858663 9.728644 



# CurvM15
options(scipen=999)
l3.thrs.CurvM15 <- kMeanThresholds(seg3.sel.sp@data$CurvM15[seg3.sel.sp@data$CurvM15 >= -1])
l3.thrs.CurvM15[order(l3.thrs.CurvM15),]
# 0.00008513638 0.00312045172 0.00671826057 0.01179443582 0.01824570166 0.02899087325 



# SfRghn
l3.thrs.SfRghn <- kMeanThresholds(seg3.sel.sp@data$SfRghn)
l3.thrs.SfRghn.3rdQu <- kMeanThresholds(seg3.sel.sp@data$SfRghn[which(seg3.sel.sp@data$SfRghn <= 0.0234500)])
l3.thrs.SfRghn[order(l3.thrs.SfRghn),]
l3.thrs.SfRghn.3rdQu[order(l3.thrs.SfRghn.3rdQu),]
# 0.001002775 0.002293221 0.003725149 0.005361888 0.007200453 0.009145321 0.011091388 0.013112721 0.015206317 0.017318034 0.019627951 0.022108847  0.033044126 0.044876192 0.059027385 0.076350798 0.102829541 0.148552598 



# SVF
l3.thrs.SVF <- kMeanThresholds(seg3.sel.sp@data$SVF)
l3.thrs.SVF[order(l3.thrs.SVF),]
# 0.7101694 0.7865066 0.8407246 0.8933890 0.9435241




# ... shape metrics --------------------------------

# P_A
l3.thrs.P_A <- kMeanThresholds(seg3.sel.sp@data$P_A)
l3.thrs.P_A[order(l3.thrs.P_A),]
# 0.2220594 0.4360765 0.7110833 1.0533773 1.4782892 3.0752165 



# D_A
l3.thrs.D_A <- kMeanThresholds(seg3.sel.sp@data$D_A)
l3.thrs.D_A[order(l3.thrs.D_A),]
# 0.03938829 0.08647260 0.14022675 0.20471096 0.28154232 0.36832040 0.47035576 0.63117312 1.13582689 



# Shp_Index
l3.thrs.Sh_Ind <- kMeanThresholds(seg3.sel.sp@data$Sh_Ind)
l3.thrs.Sh_Ind[order(l3.thrs.Sh_Ind),]
# 1.956049 2.678501 3.453671 4.472944 6.114690



# LeWiRat
l3.thrs.LeWiRat <- kMeanThresholds(seg3.sel.sp@data$LeWiRat)
l3.thrs.LeWiRat[order(l3.thrs.LeWiRat),]
# 2.523199    5.353208    8.880954   13.116816   18.194977   24.013895   30.959630   38.942666   48.166087   58.963847 73.031801   88.978868  109.061696  134.844211  168.804848  214.733364  284.377258  373.194567  485.383643  716.831079  1173.677794 4669.637236



# Convexity 
l3.thrs.Conv <- kMeanThresholds(seg3.sel.sp@data$Conv)
l3.thrs.Conv[order(l3.thrs.Conv),]
# 0.3948696 0.5888895 0.7479535 



# Compactness
l3.thrs.Comp <- kMeanThresholds(seg3.sel.sp@data$Comp)
l3.thrs.Comp[order(l3.thrs.Comp),]
# .05850803 0.11507887 0.18135528 0.26415158 0.38334628 0.61245758 





# ... textures --------------------------------

# # entropy
# flow
l3.thrs.flEntr <- kMeanThresholds(seg3.sel.sp@data$fl_Entr[which(seg3.sel.sp@data$fl_Entr != -9999)])
l3.thrs.flEntr[order(l3.thrs.flEntr),]
# 2.735676 3.545000 4.042022 4.506518 4.923196 

# object
l3.thrs.hEntr <- kMeanThresholds(seg3.sel.sp@data$h_ent_s1)
l3.thrs.hEntr[order(l3.thrs.hEntr),]
# 0.7441391 1.2101168 1.6814788 



# # IDM
# flow
l3.thrs.flIDM <- kMeanThresholds(seg3.sel.sp@data$fl_IDM[which(seg3.sel.sp@data$fl_IDM != -9999)])
l3.thrs.flIDM[order(l3.thrs.flIDM),]
# 0.1790599 0.3575048 0.5269947 

# object
l3.thrs.hIDM <- kMeanThresholds(seg3.sel.sp@data$h_idm_s1)
l3.thrs.hIDM[order(l3.thrs.hIDM),]
# 0.3566620 0.4701408 0.5600405 0.6448140 0.7234531 0.7964694 0.8872237 




# # Contrast
# flow
l3.thrs.flContr <- kMeanThresholds(seg3.sel.sp@data$fl_Contr[which(seg3.sel.sp@data$fl_Contr != -9999)])
l3.thrs.flContr.SThr <- kMeanThresholds(seg3.sel.sp@data$fl_Contr[which(seg3.sel.sp@data$fl_Contr < 19.703085 & seg3.sel.sp@data$fl_Contr != -9999)])
l3.thrs.flContr.fin <- c(l3.thrs.flContr, l3.thrs.flContr.SThr)

l3.thrs.flContr.fin[order(l3.thrs.flContr.fin)]
# 2.792844   5.282718   6.589744   8.063793  11.046769  14.340386  17.785738  19.703085  36.209643 57.049255  83.808949 117.923379 164.148080 234.003035 389.234040

# object
l3.thrs.hContr <- kMeanThresholds(seg3.sel.sp@data$h_con_s5)
l3.thrs.hContr.SThr <- kMeanThresholds(seg3.sel.sp@data$h_con_s5[which(seg3.sel.sp@data$h_con_s5 < 8.893049)])
l3.thrs.hContr.fin <- c(l3.thrs.hContr, l3.thrs.hContr.SThr)

l3.thrs.hContr.fin[order(l3.thrs.hContr.fin)]

# 0.006323301  0.436824506  0.691258307  0.975253474  1.231548163  1.500720825  1.777628495  1.846147213 2.108604472  2.526397887  2.956122668  3.410930348  3.896420398  4.228076062  4.589923137  4.933744106 5.298336167  5.680460617  6.104057967  6.520342799  6.961522945  7.455780246  8.008065372  8.575346624  8.893049 21.835414 59.492148 



# # variance
# flow
l3.thrs.flVAR <- kMeanThresholds(seg3.sel.sp@data$fl_VAR[which(seg3.sel.sp@data$fl_VAR != -9999)])
l3.thrs.flVAR.SThr <- kMeanThresholds(seg3.sel.sp@data$fl_VAR[which(seg3.sel.sp@data$fl_VAR < 27.188944 & seg3.sel.sp@data$fl_VAR != -9999)])
l3.thrs.flVAR.fin <- c(l3.thrs.flVAR[which(l3.thrs.flVAR >= 27.188944-0.1)], l3.thrs.flVAR.SThr)

l3.thrs.flVAR.fin[order(l3.thrs.flVAR.fin)]
#   1.811636   3.677751   5.754885   8.134923  10.661615  13.422412  16.277646  19.208409  22.365179  25.632585  27.188944  41.867353  58.823628  80.375781 108.990955 152.822711 239.593452 


# object
l3.thrs.hVAR <- kMeanThresholds(seg3.sel.sp@data$h_var_s1)
l3.thrs.hVAR[order(l3.thrs.hVAR),]
#  1.624923  2.707846  4.016192  5.574425  7.578808 10.135352 13.688580 19.638863 31.187728 





# ... neighbor statistics --------------------------------

# # Slope
l3.thrs.Slp.mw <- kMeanThresholds(seg3.sel.sp@data$Slp_mw[which(seg3.sel.sp@data$Slp_mw != -9999)])
l3.thrs.Slp.mw[order(l3.thrs.Slp.mw),]
# 7.496668 13.139954 19.317878 26.728982 


# mean difference in flow direction only neighbors
l3.thrs.Slp.mdfM <- kMeanThresholds(seg3.sel.sp@data$Slp_mdfM[which(seg3.sel.sp@data$Slp_mdfM != -9999)])
l3.thrs.Slp.mdfM[order(l3.thrs.Slp.mdfM),]
# 6.75608  19.09039

l3.thrs.Slp.mdfA <- kMeanThresholds(seg3.sel.sp@data$Slp_mdfA[which(seg3.sel.sp@data$Slp_mdfA != -9999)])
l3.thrs.Slp.mdfA[order(l3.thrs.Slp.mdfA),]
#  7.839448 19.541741 

l3.thrs.Slp.mwfM <- kMeanThresholds(seg3.sel.sp@data$Slp_mwfM[which(seg3.sel.sp@data$Slp_mwfM >= 0)])
l3.thrs.Slp.mwfM[order(l3.thrs.Slp.mwfM),]
# 8.157438 15.771015 24.904971 


l3.thrs.Slp.mwfA <- kMeanThresholds(seg3.sel.sp@data$Slp_mwfA[which(seg3.sel.sp@data$Slp_mwfA >= 0)])
l3.thrs.Slp.mwfA[order(l3.thrs.Slp.mwfA),]
#  9.791803 20.498408 


# # surface roughness
# mean difference in flow direction only neighbors
l3.thrs.SfRn.mdfA <- kMeanThresholds(seg3.sel.sp@data$SfRn_mdfA[which(seg3.sel.sp@data$SfRn_mdfA >= 0)])
l3.thrs.SfRn.mdfA[order(l3.thrs.SfRn.mdfA),]
# 0.002951966 0.008714301 0.014909456 0.022744592 0.033833605 0.049270337 


# weighted mean in flow direction only neighbors
l3.thrs.SfRn.mwfA <- kMeanThresholds(seg3.sel.sp@data$SfRn_mwfA[which(seg3.sel.sp@data$SfRn_mwfA >= 0)])
l3.thrs.SfRn.mwfA[order(l3.thrs.SfRn.mwfA),]
# 0.004852615 0.010766454 0.018317376 0.028887179 0.046306275 





# # openness
# mean difference in flow direction only neighbors
l3.thrs.dtmD51.mdfA <- kMeanThresholds(seg3.sel.sp@data$dtmD_mdfA[which(seg3.sel.sp@data$dtmD_mdfA != -9999)])
l3.thrs.dtmD51.mdfA[order(l3.thrs.dtmD51.mdfA),]
# mdfM: -0.3234906  0.6583425 
# -0.2529726  0.8045145 

# weighted mean in flow direction only neighbors
l3.thrs.dtmD51.mwfA <- kMeanThresholds(seg3.sel.sp@data$dtmD_mwfA[which(seg3.sel.sp@data$dtmD_mwfA != -9999)])
l3.thrs.dtmD51.mwfA[order(l3.thrs.dtmD51.mwfA),]
# mdfM: -1.0265311  0.1037361 
# -1.356766 -0.468439  0.084946 

# # entropy
# weighted mean in flow direction only neighbors
l3.thrs.fEntr.mwfM <- kMeanThresholds(seg3.sel.sp@data$fEnt_mwfM[which(seg3.sel.sp@data$fEnt_mwfM >= 0)])
l3.thrs.fEntr.mwfM[order(l3.thrs.fEntr.mwfM),]
# 3.298777 3.981127 4.543314 


l3.thrs.fEntr.mdfM  <- kMeanThresholds(seg3.sel.sp@data$fEntr_mdfM[which(seg3.sel.sp@data$fEntr_mdfM != -9999 & seg3.sel.sp@data$fEntr_mdfM  < 100)])
l3.thrs.fEntr.mdfM[order(l3.thrs.fEntr.mdfM),]
# 0.06250103 0.78817889 1.54810914 



# # Sky-view factor
# weighted mean in flow direction only neighbors
l3.thrs.fSVF.mwfM <- kMeanThresholds(seg3.sel.sp@data$SVF_mwfM[which(seg3.sel.sp@data$SVF_mwfM >= 0)])
l3.thrs.fSVF.mwfM[order(l3.thrs.fSVF.mwfM),]
# 0.8027135 0.8674388 0.9239724 

# weighted mean in flow direction only neighbors
l3.thrs.fSVF.mwfA <- kMeanThresholds(seg3.sel.sp@data$SVF_mwfA[which(seg3.sel.sp@data$SVF_mwfA >= 0)])
l3.thrs.fSVF.mwfA[order(l3.thrs.fSVF.mwfA),]
# 0.8090690 0.8742489 0.9274646 





# # contrast
# weighted mean in flow direction only neighbors
l3.thrs.fContr.mwfM <- kMeanThresholds(seg3.sel.sp@data$fCont_mwfM[which(seg3.sel.sp@data$fCont_mwfM >= 0)])
l3.thrs.fContr.mwfM[order(l3.thrs.fContr.mwfM),]
# 9.721795  19.347399  31.655446  53.312079 112.916340 

l3.thrs.fContr.mwfA <- kMeanThresholds(seg3.sel.sp@data$fCont_mwfA[which(seg3.sel.sp@data$fCont_mwfA >= 0)])
l3.thrs.fContr.mwfA[order(l3.thrs.fContr.mwfA),]
#  9.081937 17.346696 27.077867 41.364341 72.159640 


l3.thrs.fContr.mdfA <- kMeanThresholds(seg3.sel.sp@data$fCont_mdfA[which(seg3.sel.sp@data$fCont_mdfA >= -1000)])
l3.thrs.fContr.mdfA[order(l3.thrs.fContr.mdfA),]
# 38.91359  152.66443 1377.98068 2407.98566 3621.88705 5254.69941 7224.95764 9466.17864 


l3.thrs.fContr.mdfM <- kMeanThresholds(seg3.sel.sp@data$fCont_mdfM[which(seg3.sel.sp@data$fCont_mdfM != -9999 & seg3.sel.sp@data$fCont_mdfM < 81.260 )])
l3.thrs.fContr.mdfM[order(l3.thrs.fContr.mdfM),]
# -17.5465190   0.1239338  19.3817751  43.7494814  68.1475421 



# # variance
# weighted mean in flow direction only neighbors
l3.thrs.fVar.mwfM <- kMeanThresholds(seg3.sel.sp@data$fVar_mwfM[which(seg3.sel.sp@data$fVar_mwfM >= 0)])
l3.thrs.fVar.mwfM[order(l3.thrs.fVar.mwfM),]
# 8.34610 18.08403 30.52406 50.57004 92.49502 


l3.thrs.fVar.mdfA <- kMeanThresholds(seg3.sel.sp@data$fVar_mdfA[which(seg3.sel.sp@data$fVar_mdfA != -9999 & seg3.sel.sp@data$fVar_mdfA < 61.840)])
l3.thrs.fVar.mdfA[order(l3.thrs.fVar.mdfA),]
#  -12.74610   6.46147  28.79327  49.76587 


l3.thrs.fVar.mdfM <- kMeanThresholds(seg3.sel.sp@data$fVar_mdfM[which(seg3.sel.sp@data$fVar_mdfM != -9999 & seg3.sel.sp@data$fVar_mdfM < 56.550 )])
l3.thrs.fVar.mdfM[order(l3.thrs.fVar.mdfM),]
#  -19.042292  -2.210978  12.967614  30.753428  47.220238 

# -19.042292 -12.74610   -2.210978  6.46147 12.967614  30.753428  47.220238 





# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# L3: Refine Class 11 ----------------------------------------
# # # 

# ... class 41 LeWiRat intersection --------------------------
# get intersecting IDs from LeWiRat and Class 11
# subsetting and union data
L3.LWR.sub <- subset(seg3.sel.sp, seg3.sel.sp@data$LeWiRat >= 58.963847 & seg3.sel.sp@data$Class != 11 & seg3.sel.sp@data$Class != 22 & seg3.sel.sp@data$Class != 77)

# intersect geometries
L3.LWR.inters.cl11 <- unlist(gIntersects(spgeom1 = L3.LWR.sub, spgeom2 = subset(seg3.sel.sp, seg3.sel.sp@data$Class == 11), byid = TRUE, returnDense = FALSE))


# get ID's out of seg3.sel.sp
L3.classLWR.ID <- unique(subset(seg3.sel.sp, seg3.sel.sp@data$Class == 11)[L3.LWR.inters.cl11,]@data$ID)
L3.classLWR.ID.pos <- match(L3.classLWR.ID, seg3.sel.sp@data$ID)

# 554

# ... refine class 11 to 10 ------------------------------------------

# 1) "Conv" >= 0.4939739 AND "Slp_mwfM" < 24.904971 OR "Conv" >= 0.4939739  AND "DTM_NH" < 0.53951430
L3.classLWR.Q <-  Reduce(intersect, list(L3.classLWR.ID.pos, which(seg3.sel.sp@data$Class == 11 & seg3.sel.sp@data$Conv >= 0.4939739 & (seg3.sel.sp@data$Slp_mwfM < 24.904971 |  seg3.sel.sp@data$DTM_NH < 0.53951430))))

# set class 11 to 10
seg3.sel.sp@data$Class[which(seg3.sel.sp@data$Class == 11)] <- 10

# kick out L3.classLWR.Q
seg3.sel.sp@data$Class[L3.classLWR.Q] <- 11

# 2) "Slp_mdfM" < 0 OR "Slp_mdfA"  < 0 OR "Slp_mdfM" < 6.756081 AND "Conv" >= 0.4975 AND "Slp_mwfM" >= 0  AND  "fl_Entr" < 4.923196 - mean difference slope
# 3) "SVF_mwfM" >= 0.9239724 AND "Conv" >= 0.497 OR "SVF_mwfA" >= 0.9239724 AND "Slp_mwfA" < 6.756081 - sky view factor
# 4) "NH_rB_11_s" > 40 or "NH_rB_11_l"  > 75 -  relative border to NH
# 5) "fVar_mdfA" >= 60  AND "Slp_mwfA" <= 9.791803 - Step 5: Earth and river banks
# 6) "P_A" >= 0.7110833 AND "Comp" <= 0.14756019   AND  "fVar_mdfA" >= 80 - Step 6: Again small elements
# 7) "Area" >= 850 AND "fEntr_mdfM" < 0 AND "Slp_mw" >= 13.139954 OR "LeWiRat" >= 18.194977  AND "fEntr_mdfM" < 0 OR "Area" >= 350 AND "NH_rB_11_s" > 0 AND "NH_rB_11_l" > 25 - large areas
# 8) "LWR_rB_11" > 10  AND "NH_rB_11_l" > 0 - neighborhood to LWR + NH
seg3.sel.sp@data$Class[which(seg3.sel.sp@data$Class == 10 & ((seg3.sel.sp@data$Slp_mdfM < 0 | seg3.sel.sp@data$Slp_mdfA < 0 | seg3.sel.sp@data$Slp_mdfM < 6.756081 & seg3.sel.sp@data$Conv >= 0.4939739 & seg3.sel.sp@data$Slp_mwfM >= 0 & seg3.sel.sp@data$fl_Entr < 4.923196)
                             | (seg3.sel.sp@data$SVF_mwfM >= 0.9239724 & seg3.sel.sp@data$Conv >= 0.4939739 | seg3.sel.sp@data$SVF_mwfA >= 0.9274646 & seg3.sel.sp@data$Slp_mwfA < 6.756081)
                             | (seg3.sel.sp@data$NH_rB_11_s >= 40 | seg3.sel.sp@data$NH_rB_11_l >= 75 | seg3.sel.sp@data$LWR_rB_11 >= 10 & seg3.sel.sp@data$NH_rB_11_l > 0)
                             | (seg3.sel.sp@data$fVar_mdfA >= 60 & seg3.sel.sp@data$Slp_mwfA <= 9.791803)
                             | (seg3.sel.sp@data$P_A >= 0.7110833 & seg3.sel.sp@data$Comp < 0.14756019 & seg3.sel.sp@data$fVar_mdfA >= 80)
                             | (seg3.sel.sp@data$Area >= 850 & seg3.sel.sp@data$fEntr_mdfM < 0 & seg3.sel.sp@data$Slp_mw >= 13.139954 | seg3.sel.sp@data$LeWiRat >= 18.194977 & seg3.sel.sp@data$fEntr_mdfM < 0 | seg3.sel.sp@data$Area >= 350 & seg3.sel.sp@data$NH_rB_11_s > 0 & seg3.sel.sp@data$NH_rB_11_l >= 25)))] <- 11

# 1736 objects could be removed

# length(which(seg3.sel.sp@data$Class == 10 & (seg3.sel.sp@data$Slp_mdfM < 0 | (seg3.sel.sp@data$Slp_mdfA < 0 | seg3.sel.sp@data$Slp_mdfM < 6.756081 & seg3.sel.sp@data$Conv >= 0.4939739 & seg3.sel.sp@data$Slp_mwfM >= 0 & seg3.sel.sp@data$fl_Entr < 4.923196)))) # 512
# length(which(seg3.sel.sp@data$Class == 10 & ((seg3.sel.sp@data$SVF_mwfM >= 0.9239724 & seg3.sel.sp@data$Conv >= 0.4939739) | (seg3.sel.sp@data$SVF_mwfA >= 0.9239724 & seg3.sel.sp@data$Slp_mwfA < 6.756081))))
# length(which(seg3.sel.sp@data$Class == 10 & (seg3.sel.sp@data$NH_rB_11_s >= 40 | seg3.sel.sp@data$NH_rB_11_l >= 75 | seg3.sel.sp@data$LWR_rB_11 >= 10 & seg3.sel.sp@data$NH_rB_11_l > 0)))
# length(which(seg3.sel.sp@data$Class == 10 & (seg3.sel.sp@data$fVar_mdfA >= 60 & seg3.sel.sp@data$Slp_mwfA <= 9.791803)))
# length(which(seg3.sel.sp@data$Class == 10 & (seg3.sel.sp@data$P_A >= 0.7110833 & seg3.sel.sp@data$Comp < 0.14756019 & seg3.sel.sp@data$fVar_mdfA >= 80)))
# length(which(seg3.sel.sp@data$Class == 10 & ((seg3.sel.sp@data$Area >= 850 & seg3.sel.sp@data$fEntr_mdfM < 0 & seg3.sel.sp@data$Slp_mw >= 13.139954) | (seg3.sel.sp@data$LeWiRat >= 18.194977 & seg3.sel.sp@data$fEntr_mdfM < 0) | (seg3.sel.sp@data$Area >= 350 & seg3.sel.sp@data$NH_rB_11_s > 0 & seg3.sel.sp@data$NH_rB_11_s >= 25))))



# ... refine class 22 and 77 as well ----------------------------
# intersection with class 10
L3.class10 <- subset(seg3.sel.sp, seg3.sel.sp@data$Class == 10)

# intersect geometries
L3.cl22.inters.cl10 <- unlist(gIntersects(spgeom1 = L3.class10, spgeom2 = subset(seg3.sel.sp, seg3.sel.sp@data$Class == 22), byid = TRUE, returnDense = FALSE))
L3.cl77.inters.cl10 <- unlist(gIntersects(spgeom1 = L3.class10, spgeom2 = subset(seg3.sel.sp, seg3.sel.sp@data$Class == 77), byid = TRUE, returnDense = FALSE))



# get ID's out of L1.class111 subset and get position afterwards
L3.cl22.inters.cl10.ID <- subset(seg3.sel.sp, seg3.sel.sp@data$Class == 22)[L3.cl22.inters.cl10,]@data$ID
L3.cl22.inters.cl10.pos <- match(L3.cl22.inters.cl10.ID, seg3.sel.sp@data$ID)

L3.cl77.inters.cl10.ID <- subset(seg3.sel.sp, seg3.sel.sp@data$Class == 77)[L3.cl77.inters.cl10,]@data$ID
L3.cl77.inters.cl10.pos <- match(L3.cl77.inters.cl10.ID, seg3.sel.sp@data$ID)


# set class 20 and 70 into class 10
seg3.sel.sp@data$Class[L3.cl22.inters.cl10.pos] <- 20
seg3.sel.sp@data$Class[L3.cl77.inters.cl10.pos] <- 70


# get class objects in L3
L3.class10.pos <- which(seg3.sel.sp@data$Class == 10)
L3.class20.pos <- which(seg3.sel.sp@data$Class == 20)
L3.class70.pos <- which(seg3.sel.sp@data$Class == 70)


# get neighbors of class objects in L3
L3.nb.class10 <- unique(unlist(L3.nb[L3.class10.pos]))
L3.nb.class20 <- unique(unlist(L3.nb[L3.class20.pos]))
L3.nb.class70 <- unique(unlist(L3.nb[L3.class70.pos]))







# # #
# L3: Landslide Body for Class 10--------------------------------------------
# # #

# ... refine bounding box and neighbors in flow direction ----------------------------

# ID name is identical to subset(seg3.sel.sp, seg3.sel.sp@data$Class == 11)@data$ID
L3.bb.flow.cl10.pos <- match(seg3.sel.sp@data$ID[which(seg3.sel.sp@data$Class == 10)], L3.bb.flow.cl11$ID)
L3.bb.flow.cl10  <- L3.bb.flow.cl11[L3.bb.flow.cl10.pos,]

writeSpatialShape(x = L3.bb.flow.cl10, fn = paste(getwd(), dirname(seg3.sel), "L3_bb_flow_cl10.shp", sep = "/"))



# ... ... get neighbors in flow direction ---------------------------------
# get all objects that intersect the bounding boxes
L3.bbFlowCl10.inters.seg3 <- unique(unlist(gIntersects(spgeom1 = L3.bb.flow.cl10, spgeom2 = seg3.sel.sp, byid = TRUE, returnDense = FALSE)))


# kick out class objects
L3.bbFlowCl10.inters.seg3 <- L3.bbFlowCl10.inters.seg3[!L3.bbFlowCl10.inters.seg3 %in% c(L3.class10.pos, L3.class20.pos, L3.class70.pos)] # position of L3.class11.pos has not changed


# kick out not neighbor objects
L3.bbFlowCl10.inters.seg3.nb  <- L3.bbFlowCl10.inters.seg3[L3.bbFlowCl10.inters.seg3 %in% c(L3.nb.class10)]


# put information to data 
seg3.sel.sp@data$Cl10Fl <- 0
seg3.sel.sp@data$Cl10FlNb <- 0

seg3.sel.sp$Cl10Fl[L3.bbFlowCl10.inters.seg3] <- 1
seg3.sel.sp$Cl10FlNb[L3.bbFlowCl10.inters.seg3.nb] <- 1



# ... RelationalClass Function ------------------------
# statistics of neighbors to a given class

# surface roughness
# [1] "------ Run of RelationalClassFunction: 12.4676666666667 Minutes ------"
L3.SfRghn.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("SfRghn", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.SfRghn.RelClass.10 <- L3.SfRghn.RelClass[c("ID", "m_w", "m_d_nb_cl")]
colnames(L3.SfRghn.RelClass.10) <- c("ID", "SfRn_mw2Cl", "SfRn_md2Cl")




# DTM
L3.DTM.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("DTM", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.DTM.RelClass.10 <- L3.DTM.RelClass [c("ID", "bor_cl_rel", "m_w", "m_d_nb_cl")]
colnames(L3.DTM.RelClass.10) <- c("ID", "rBor_Cl10", "DTM_mw2Cl", "DTM_md2Cl")




# sky-view factor
L3.SVF.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("SVF", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.SVF.RelClass.10 <- L3.SVF.RelClass [c("ID", "m_w", "m_d_nb_cl")]
colnames(L3.SVF.RelClass.10) <- c("ID", "SVF_mw2Cl", "SVF_md2Cl")


# slope
L3.Slp.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("Slp", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.Slp.RelClass.10 <- L3.Slp.RelClass [c("ID", "m_w", "m_d_nb_cl")]
colnames(L3.Slp.RelClass.10) <- c("ID", "Slp_mw2Cl", "Slp_md2Cl")




# slope SD
L3.SlpSD.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("Slp_SD", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.SlpSD.RelClass.10 <- L3.SlpSD.RelClass [c("ID", "m_w", "m_d_nb_cl")]
colnames(L3.SlpSD.RelClass.10) <- c("ID", "SlpSD_mw2Cl", "SlpSD_md2Cl")



# opennness
L3.dtmD51.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("DTM_D_51", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.dtmD51.RelClass.10 <- L3.dtmD51.RelClass [c("ID", "m_w", "m_d_nb_cl")]
colnames(L3.dtmD51.RelClass.10) <- c("ID", "dtmD_mw2Cl", "dtmD_md2Cl")



# entropy flow
L3.flEntr.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("fl_Entr", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.flEntr.RelClass.10 <- L3.flEntr.RelClass [c("ID", "m_w", "m_d_nb_cl")]
colnames(L3.flEntr.RelClass.10) <- c("ID", "fEnt_mw2Cl", "fEnt_md2Cl")



# entropy object
L3.hEntr.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("h_ent_s1", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.hEntr.RelClass.10 <- L3.hEntr.RelClass [c("ID",  "m_w", "m_d_nb_cl")]
colnames(L3.hEntr.RelClass.10) <- c("ID",  "hEnt_mw2Cl", "hEnt_md2Cl")



# IDM flow
L3.flIDM.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("fl_IDM", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.flIDM.RelClass.10 <- L3.flIDM.RelClass [c("ID", "m_w", "m_d_nb_cl")]
colnames(L3.flIDM.RelClass.10) <- c("ID",  "fIDM_mw2Cl", "fIDM_md2Cl")



# contrast flow
L3.flContr.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("fl_Contr", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.flContr.RelClass.10 <- L3.flContr.RelClass [c("ID",  "m_w", "m_d_nb_cl")]
colnames(L3.flContr.RelClass.10) <- c("ID",  "fContr_mw2Cl", "fContr_md2Cl")




# contrast object
L3.hContr.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("h_con_s5", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.hContr.RelClass.10 <- L3.hContr.RelClass [c("ID",  "m_w", "m_d_nb_cl")]
colnames(L3.hContr.RelClass.10) <- c("ID",  "hContr_mw2Cl", "hContr_md2Cl")



# variance flow
L3.flVar.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("fl_VAR", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.flVar.RelClass.10 <- L3.flVar.RelClass [c("ID",  "m_w", "m_d_nb_cl")]
colnames(L3.flVar.RelClass.10) <- c("ID",  "fVar_mw2Cl", "fVar_md2Cl")




# LeWiRat
L3.LWR.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("LeWiRat", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.LWR.RelClass.10 <- L3.LWR.RelClass [c("ID",  "m_w", "m_d_nb_cl", "rat_nb_cl")]
colnames(L3.LWR.RelClass.10) <- c("ID",  "LWR_mw2Cl", "LWR_md2Cl", "LWR_r2Cl")




# Area
L3.Area.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("Area", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.Area.RelClass.10 <- L3.Area.RelClass [c("ID", "m_d_nb_cl", "rat_nb_cl")]
colnames(L3.Area.RelClass.10) <- c("ID", "A_md2Cl", "A_r2Cl")



# Shape Index
L3.ShInd.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb, class.var = 10, var = c("Sh_Ind", "Area"), quiet = FALSE, class.as.neighbors = FALSE, calc.nb.flow = TRUE, bb = L3.bb.flow.cl10)
L3.ShInd.RelClass.10 <- L3.ShInd.RelClass [c("ID", "m_w", "m_d_nb_cl", "rat_nb_cl")]
colnames(L3.ShInd.RelClass.10) <- c("ID", "SId_mw2Cl", "SId_md2Cl", "SId_r2Cl")




# ...  merge to table -----------------------------------
L3.RelClass.10 <- cbind(L3.DTM.RelClass.10, L3.Slp.RelClass.10[-1], L3.SVF.RelClass.10[-1], L3.flVar.RelClass.10[-1], L3.flEntr.RelClass.10[-1], L3.LWR.RelClass.10[-1], L3.Area.RelClass.10[-1], L3.ShInd.RelClass.10[-1])
seg3.sel.sp@data <- merge(seg3.sel.sp@data, L3.RelClass.10 , by = "ID", all.x = TRUE)



# seg3.sel.sp@data <- seg3.sel.sp@data[-c( (ncol(seg3.sel.sp@data) - 17): ncol(seg3.sel.sp@data))]




# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ... calculate k-mean thresholds ----------------------------------------
# # # 

# # Slope
l3.thrs.Slp.md2Cl <- kMeanThresholds(seg3.sel.sp@data$Slp_md2Cl[which(seg3.sel.sp@data$Slp_md2Cl >= 0)])
l3.thrs.Slp.md2Cl [order(l3.thrs.Slp.md2Cl),]
#  1.409017  4.336116  7.758105 12.848418 


l3.thrs.Slp.md2Cl.N <- kMeanThresholds(seg3.sel.sp@data$Slp_md2Cl[which(seg3.sel.sp@data$Slp_md2Cl < 0 & seg3.sel.sp@data$Slp_md2Cl != -9999)])
l3.thrs.Slp.md2Cl.N [order(l3.thrs.Slp.md2Cl.N),]
# -21.073871 -13.309070  -5.491192 


# # Slope
l3.thrs.Slp.mw2Cl <- kMeanThresholds(seg3.sel.sp@data$Slp_mw2Cl[which(seg3.sel.sp@data$Slp_mw2Cl >= 0)], size.sample = NULL)
l3.thrs.Slp.mw2Cl[order(l3.thrs.Slp.mw2Cl),]
#  14.04851 26.54378 



# # DTM
l3.thrs.DTM.md2Cl <- kMeanThresholds(seg3.sel.sp@data$DTM_md2Cl[which(seg3.sel.sp@data$DTM_md2Cl > -9999 & seg3.sel.sp@data$DTM_md2Cl <= 0)])
l3.thrs.DTM.md2Cl[order(l3.thrs.DTM.md2Cl),]
#  15.476275 -10.767940  -7.554506  -4.931933  -2.826722  -0.937815  0.916749  2.840925  4.913497  7.415527 10.936923 17.214348





# # Shape Index
l3.thrs.ShInd.mw2Cl <- kMeanThresholds(seg3.sel.sp@data$SId_mw2Cl[which(seg3.sel.sp@data$SId_mw2Cl >= 0)])
l3.thrs.ShInd.mw2Cl[order(l3.thrs.ShInd.mw2Cl),]
#  2.601358 3.705090 5.228935 



# # Length Width Ratio
l3.thrs.LWR.mw2Cl <- kMeanThresholds(seg3.sel.sp@data$LWR_mw2Cl[which(seg3.sel.sp@data$LWR_mw2Cl >= 0)])
l3.thrs.LWR.mw2Cl[order(l3.thrs.LWR.mw2Cl),]
#  5.145388  12.971049  24.186127  43.664810  82.280344 148.080811 463.008127 



# # Area
l3.thrs.Area.r2Cl <- kMeanThresholds(seg3.sel.sp@data$A_r2Cl[which(seg3.sel.sp@data$A_r2Cl >= 0)])
l3.thrs.Area.r2Cl[order(l3.thrs.Area.r2Cl),]
#  1.131998   6.859520  14.473987  24.306542  37.315021  54.347445  78.131136 115.310797 166.458431 253.175149 392.809488 

l3.thrs.Area.r2Cl.S <- kMeanThresholds(seg3.sel.sp@data$A_r2Cl[which(seg3.sel.sp@data$A_r2Cl >= 0 & seg3.sel.sp@data$A_r2Cl < 12.0100)])
l3.thrs.Area.r2Cl.S[order(l3.thrs.Area.r2Cl.S),]
#  0.1773351  0.7331048  1.3596988  2.0977633  2.9970024  3.9494035  5.1864926  6.5547397  8.0506983  9.4997718 11.1052252



# # Border
l3.thrs.relBor <- kMeanThresholds(seg3.sel.sp@data$rBor_Cl10[which(seg3.sel.sp@data$rBor_Cl10 >= 0)])
l3.thrs.relBor[order(l3.thrs.relBor),]
#  5.977698 21.815241 44.763842 93.617249 




# 1.) "fl_Entr" >= 4.042022 AND "SVF" < 0.9435241 --> 4854
# 2.) "Slp_md2Cl" < 1.409017 AND "DTM_md2Cl" < 2.840925 AND "SfRghn" >= 0.003725149 
# 3.) "Cl10FlNb" = 0 AND "Area" >= 2000 OR "Cl10FlNb" = 0 AND "LeWiRat" >= 24.013895 OR "Cl10FlNb" = 0 AND "Sh_Ind"  >= 4.47294 OR "Cl10FlNb" = 0 AND "Slp" >=  27.504224 OR "Cl10FlNb" = 0 AND "DTM_D_51" <   -3.2177916 OR "Cl10FlNb" = 0 AND "DTM_NH" < 0.08239538
# 4.) "Slp_mw2Cl"  >=  26.54378 AND "Area" >= 150 AND "rBor_Cl10"  < 10 OR  "Slp_md2Cl" > 0 AND "Slp_mw2Cl" >= 26.54378
# 5.) "A_r2Cl" >= 6.5547397  AND  ("Slp" >= 27.504224 OR "Slp" < 8.821035  OR "SfRghn" < 0.005361888 OR  "DTM_NH" < 0.19024493  AND "rBor_Cl10" < 5.977698 OR "DTM_md2Cl" > 0 OR "Slp_md2Cl" > 0)
# 6.) "A_md2Cl" < 0 AND "A_md2Cl" > -9999 AND "A_r2Cl" < 0.7331048  AND "rBor_Cl10" < 15 OR "A_r2Cl" < 0.1 AND  "A_r2Cl" > 0 AND "rBor_Cl10" < 50 OR  "A_r2Cl" < 1 AND "A_md2Cl" > -9999 AND ("Slp_md2Cl" > 0 OR "DTM_md2Cl" > 0)
# 7.) "SId_mw2Cl"  >= 3.705090 AND "rBor_Cl10" < 5 OR "SId_mw2Cl"  >= 5.228935 OR "SId_mw2Cl"  >= 3.705090 AND "DTM_NH" < 0.08239538 OR "SId_mw2Cl"  >= 3.705090 AND "Slp" < 8.821035 OR "SId_mw2Cl"  >= 3.705090 AND "Slp"> 26.54378  AND "DTM_md2Cl" < -0.937815  OR "SId_mw2Cl"  >=  4.47294 AND "A_r2Cl" >= 5.1864926
# 8.) "LWR_mw2Cl" >= 12.971049 AND "Slp" <= 8.821035  AND ("rBor_Cl10" <= 29 OR "SId_md2Cl" > 0 OR "DTM_md2Cl" > 0) OR "LWR_mw2Cl" >=  30.959630

# ... GET potential Landslide Body ------------------------
seg3.sel.sp$Class[which(seg3.sel.sp@data$Class != 10 & seg3.sel.sp@data$Class != 11 & seg3.sel.sp@data$Class != 20 & seg3.sel.sp@data$Class != 22 & seg3.sel.sp@data$Class != 70 & seg3.sel.sp@data$Class != 77 & seg3.sel.sp@data$Cl10Fl == 1
                        & ((seg3.sel.sp@data$fl_Entr >= 4.042022 & seg3.sel.sp@data$SVF < 0.9435241)
                        & (seg3.sel.sp@data$Slp_md2Cl < 1.409017 & seg3.sel.sp@data$DTM_md2Cl < 2.840925 & seg3.sel.sp@data$SfRghn >= 0.003725149)
                        & !(seg3.sel.sp@data$Cl10FlNb == 0 & (seg3.sel.sp@data$Area >= 2000 | seg3.sel.sp@data$LeWiRat >= 24.013895 | seg3.sel.sp@data$Sh_Ind  >= 4.47294 | seg3.sel.sp@data$Slp >=  27.504224 | seg3.sel.sp@data$DTM_D_51 <   -3.2177916 | seg3.sel.sp@data$DTM_NH < 0.08239538))
                        & !(seg3.sel.sp@data$Slp_mw2Cl >=  26.54378 & seg3.sel.sp@data$Area >= 150 & seg3.sel.sp@data$rBor_Cl10 < 10 | seg3.sel.sp@data$Slp_md2Cl > 0 & seg3.sel.sp@data$Slp_mw2Cl >= 26.54378)
                        & !(seg3.sel.sp@data$A_r2Cl >= 6.5547397 & (seg3.sel.sp@data$Slp >= 27.504224 | seg3.sel.sp@data$Slp < 8.821035 | seg3.sel.sp@data$SfRghn < 0.005361888 | seg3.sel.sp@data$DTM_NH < 0.19024493 & seg3.sel.sp@data$rBor_Cl10 < 5.977698 | seg3.sel.sp@data$DTM_md2Cl > 0 | seg3.sel.sp@data$Slp_md2Cl > 0) | seg3.sel.sp@data$A_r2Cl >= 37.315021 | seg3.sel.sp@data$A_md2Cl >= 5000)
                        & !(seg3.sel.sp@data$A_md2Cl > -9999 & ((seg3.sel.sp@data$A_md2Cl < 0 & seg3.sel.sp@data$A_r2Cl < 0.7331048 & seg3.sel.sp@data$rBor_Cl10 < 15) | (seg3.sel.sp@data$A_r2Cl < 0.1 & seg3.sel.sp@data$rBor_Cl10 < 50)  | (seg3.sel.sp@data$A_r2Cl < 1 & (seg3.sel.sp@data$Slp_md2Cl > 0 | seg3.sel.sp@data$DTM_md2Cl > 0))))
                        & !(seg3.sel.sp@data$SId_mw2Cl >= 3.705090 & (seg3.sel.sp@data$rBor_Cl10 < 5 | seg3.sel.sp@data$DTM_NH < 0.08239538 | seg3.sel.sp@data$Slp < 8.821035 | seg3.sel.sp@data$Slp >= 26.54378 & seg3.sel.sp@data$DTM_md2Cl < -0.937815 ) | seg3.sel.sp@data$SId_mw2Cl >= 5.22893 | seg3.sel.sp@data$SId_mw2Cl >= 4.47294 & seg3.sel.sp@data$A_r2Cl >= 5.1864926)
                        & !(seg3.sel.sp@data$LWR_mw2Cl >= 12.971049 & seg3.sel.sp@data$Slp < 8.821035 & (seg3.sel.sp@data$rBor_Cl10 < 29 | seg3.sel.sp@data$DTM_md2Cl > 0 | seg3.sel.sp@data$Slp_md2Cl > 0) | seg3.sel.sp@data$LWR_mw2Cl >=  30.959630)
                        & !(seg3.sel.sp@data$rBor_Cl10 == 100)
                        ))] <- 30

# 1792
# length(seg3.sel.sp$Class[which(seg3.sel.sp@data$Class != 10 & seg3.sel.sp@data$Class != 11 & seg3.sel.sp@data$Class != 20 & seg3.sel.sp@data$Class != 22 & seg3.sel.sp@data$Class != 70 & seg3.sel.sp@data$Class != 77 & seg3.sel.sp@data$Cl10Fl == 1)])
# 6503

# ... Class 1 ---------------------------------
L3.class30 <- subset(seg3.sel.sp, seg3.sel.sp@data$Class ==  30)


L3.class30.inters.class10 <- unique(unlist(gIntersects(spgeom1 = L3.class30, spgeom2 = L3.class10, byid = TRUE, returnDense = FALSE)))

# get ID's out of L3.class10 subset and get position afterwards
L3.class30.inters.class10.ID <- L3.class10[L3.class30.inters.class10,]@data$ID
L3.class30.inters.class10.pos <- match(L3.class30.inters.class10.ID , seg3.sel.sp@data$ID)

# set class 1
seg3.sel.sp@data$Class[L3.class30.inters.class10.pos] <- 1



# ... neighbors in inverse flow direction ---------------------------------
# calculation of bounding box in flow and inverse flow direction
L3.bb.inv.flow.Conv <- getBoundingBox(shape = subset(seg3.sel.sp, seg3.sel.sp@data$Conv < 0.4667959 & seg3.sel.sp@data$Class == 1), scale.factor = c(0.3, 0.3), 
                                  k.centroid = 2, set.centroid = "direction", scale.side = "small", centroid = FALSE, col.name = "FlowInv", quiet = TRUE)

row.names(L3.bb.inv.flow.Conv) <- row.names(subset(seg3.sel.sp, seg3.sel.sp@data$Conv < 0.4667959 & seg3.sel.sp@data$Class == 1))


L3.bb.inv.flow.NotConv <- getBoundingBox(shape = subset(seg3.sel.sp, seg3.sel.sp@data$Conv >= 0.4667959 & seg3.sel.sp@data$Class == 1), scale.factor = c(0.3, 0.3), 
                                      k.centroid = 4, set.centroid = "direction", scale.side = "small", centroid = FALSE, col.name = "FlowInv", quiet = TRUE)

row.names(L3.bb.inv.flow.NotConv) <- row.names(subset(seg3.sel.sp, seg3.sel.sp@data$Conv >= 0.4667959 & seg3.sel.sp@data$Class == 1))


# merge both bounding boxes
L3.bb.flow.inv.cl1 <- rbind(L3.bb.inv.flow.Conv, L3.bb.inv.flow.NotConv)

# order bounding box to fit ID classes order
L3.bb.flow.inv.cl1  <- L3.bb.flow.inv.cl1[order(L3.bb.flow.inv.cl1$ID),]



writeSpatialShape(x = L3.bb.flow.inv.cl1, fn = paste(getwd(), dirname(seg3.sel), "L3_bb_flow_inv_cl1.shp", sep = "/"))


# get all objects that intersect the bounding boxes
L3.bbFlowInvCl1.inters.seg3 <- unique(unlist(gIntersects(spgeom1 = L3.bb.flow.inv.cl1, spgeom2 = seg3.sel.sp, byid = TRUE, returnDense = FALSE)))


# kick out class objects
L3.bbFlowInvCl1.inters.seg3 <- L3.bbFlowInvCl1.inters.seg3[!L3.bbFlowInvCl1.inters.seg3 %in% c(L3.class10.pos, L3.class20.pos, L3.class70.pos)]


# put information to data 
seg3.sel.sp$Cl1FlInv <- 0

seg3.sel.sp$Cl1FlInv[L3.bbFlowInvCl1.inters.seg3] <- 1






# ... ClassNeighbor Function --------------------------------------------------
# Class to Class To expression
# [1] "------ Run of ClassNeighborFunction: 1.10516666666666 Minutes ------"

# border of normalized height to class 11
# 0.08239538 0.19024493
# ... ... ... border between Class 30 and Class 1 ---------------------------------
L3.Cl30.RelClass.Cl1 <- RelationalClassFunction(spdf = seg3.sel.sp, nb =  L3.nb,  class.var = 30, class.var2 = 1, var = c("Slp", "Area"), quiet = FALSE, class.as.neighbors = FALSE)

L3.Cl30.RelClass.Cl1 <- L3.Cl30.RelClass.Cl1[which(L3.Cl30.RelClass.Cl1$Cl == 1),][c("ID", "bor_cl_rel")]

L3.Cl30.RelClass.Cl1[is.na(L3.Cl30.RelClass.Cl1)] <- 0
colnames(L3.Cl30.RelClass.Cl1) <- c("ID", "rB_Cl_30_1")








# ...  merge to table -----------------------------------
seg3.sel.sp@data <- merge(seg3.sel.sp@data, L3.Cl30.RelClass.Cl1 , by = "ID", all.x = TRUE)


# seg3.sel.sp@data <- seg3.sel.sp@data[-c( (ncol(seg3.sel.sp@data) - 17): ncol(seg3.sel.sp@data))]





# ... Last Refinements before Class 3  --------------------------------------------------
# Class 1
# "rB_Cl_30_1" = 100 OR "rB_Cl_30_1" >=  90 AND "Conv" >= 0.4939739 OR "rB_Cl_30_1" = 0 OR "rB_Cl_30_1" < 14 AND "Sh_Ind" >= 3.8
seg3.sel.sp$Class[which(seg3.sel.sp@data$Class == 1 & (seg3.sel.sp@data$rB_Cl_30_1 == 100 | seg3.sel.sp@data$rB_Cl_30_1 == 0
                                                       | (seg3.sel.sp@data$rB_Cl_30_1 >= 90 & seg3.sel.sp@data$Conv >= 0.4939739)
                                                       | (seg3.sel.sp@data$rB_Cl_30_1 < 14 & seg3.sel.sp@data$Sh_Ind >= 3.8)))] <- 10


# Class 30 
# "Cl1FlInv" = 1 AND  "A_r2Cl" >= 24.306542 OR "Cl1FlInv" = 1 AND "DTM_md2Cl" > 0 AND "rBor_Cl10" < 30
seg3.sel.sp$Class[which(seg3.sel.sp@data$Class == 30 & seg3.sel.sp@data$Cl1FlInv == 1
                        & (seg3.sel.sp@data$A_r2Cl >= 24.306542 | seg3.sel.sp@data$DTM_md2Cl > 0 & seg3.sel.sp@data$rBor_Cl10 < 30))] <- -9999



# recheck neighbors again
L3.class1 <- subset(seg3.sel.sp, seg3.sel.sp@data$Class == 1)
L3.class30 <- subset(seg3.sel.sp, seg3.sel.sp@data$Class ==  30) # overwrite
L3.class30.pos <- which(seg3.sel.sp@data$Class == 30) 

L3.class30.inters.class1 <- unique(unlist(gIntersects(spgeom1 = L3.class30, spgeom2 = L3.class1, byid = TRUE, returnDense = FALSE)))

# get ID's out of L1.class111 subset and get position afterwards
L3.class30.inters.class1.ID <- L3.class1[L3.class30.inters.class1,]@data$ID
L3.class30.inters.class1.pos <- match(L3.class30.inters.class1.ID , seg3.sel.sp@data$ID)

# set class 10 and 1
seg3.sel.sp@data$Class[which(seg3.sel.sp@data$Class == 1)] <- 10
seg3.sel.sp@data$Class[L3.class30.inters.class1.pos] <- 1






# # #
# L3: Landslide Body - Class 3 --------------------------------------------
# # #



# ... Get Flanks and Crowns of Scarps (Class 2 and 7) --------------------------------------------------

# overwrite class 1 data
L3.class1 <- subset(seg3.sel.sp, seg3.sel.sp@data$Class == 1)

# intersect geometries
L3.cl20.inters.cl1 <- unlist(gIntersects(spgeom1 = L3.class1, spgeom2 = subset(seg3.sel.sp, seg3.sel.sp@data$Class == 20), byid = TRUE, returnDense = FALSE))
L3.cl70.inters.cl1 <- unlist(gIntersects(spgeom1 = L3.class1, spgeom2 = subset(seg3.sel.sp, seg3.sel.sp@data$Class == 70), byid = TRUE, returnDense = FALSE))


# get ID's out of L1.class1 subset and get position afterwards
L3.cl20.inters.cl1.ID <- subset(seg3.sel.sp, seg3.sel.sp@data$Class == 20)[L3.cl20.inters.cl1,]@data$ID
L3.cl20.inters.cl1.pos <- match(L3.cl20.inters.cl1.ID, seg3.sel.sp@data$ID)

L3.cl70.inters.cl1.ID <- subset(seg3.sel.sp, seg3.sel.sp@data$Class == 70)[L3.cl70.inters.cl1,]@data$ID
L3.cl70.inters.cl1.pos <- match(L3.cl70.inters.cl1.ID, seg3.sel.sp@data$ID)


# set class 2 and 7
seg3.sel.sp@data$Class[L3.cl20.inters.cl1.pos] <- 2
seg3.sel.sp@data$Class[L3.cl70.inters.cl1.pos] <- 7


# get class objects in L3
L3.class1.pos <- which(seg3.sel.sp@data$Class == 1)
L3.class2.pos <- which(seg3.sel.sp@data$Class == 2) # 116
L3.class7.pos <- which(seg3.sel.sp@data$Class == 7) # 304


# get neighbors of class objects in L3
L3.nb.class1 <- unique(unlist(L3.nb[L3.class1.pos]))
L3.nb.class2 <- unique(unlist(L3.nb[L3.class2.pos]))
L3.nb.class7 <- unique(unlist(L3.nb[L3.class7.pos]))



# ... Class 1: refine bounding box and neighbors in flow direction ----------------------------

# ID name is identical to subset(seg3.sel.sp, seg3.sel.sp@data$Class == 11)@data$ID
L3.bb.flow.cl1.pos <- match(seg3.sel.sp@data$ID[which(seg3.sel.sp@data$Class == 1)], L3.bb.flow.cl11$ID)
L3.bb.flow.cl1  <- L3.bb.flow.cl11[L3.bb.flow.cl1.pos,]

writeSpatialShape(x = L3.bb.flow.cl1, fn = paste(getwd(), dirname(seg3.sel), "L3_bb_flow_cl1.shp", sep = "/"))



# ... ... get neighbors in flow direction ---------------------------------
# get all objects that intersect the bounding boxes
L3.bbFlowCl1.inters.seg3 <- unique(unlist(gIntersects(spgeom1 = L3.bb.flow.cl1, spgeom2 = seg3.sel.sp, byid = TRUE, returnDense = FALSE)))


# kick out class objects
L3.bbFlowCl1.inters.seg3 <- L3.bbFlowCl1.inters.seg3[!L3.bbFlowCl1.inters.seg3 %in% c(L3.class1.pos, L3.class2.pos, L3.class7.pos, L3.class30.pos)] # position of L3.class11.pos has not changed


# kick out not neighbor objects
L3.bbFlowCl1.inters.seg3.nb  <- L3.bbFlowCl1.inters.seg3[L3.bbFlowCl1.inters.seg3 %in% c(L3.nb.class1, L3.nb.class2)]


# put information to data 
seg3.sel.sp@data$Cl1Fl <- 0
seg3.sel.sp@data$Cl1FlNb <- 0

seg3.sel.sp$Cl1Fl[L3.bbFlowCl1.inters.seg3] <- 1
seg3.sel.sp$Cl1FlNb[L3.bbFlowCl1.inters.seg3.nb] <- 1





# ... border of flow polygons to Class 1 and 30 ----------------------------

# create data for analyses
seg3.sel.sp.body <- seg3.sel.sp

seg3.sel.sp.body@data$Class[which(seg3.sel.sp@data$Class == 1 | seg3.sel.sp@data$Class == 30)] <- 130 # 2376
seg3.sel.sp.body@data$Class[which(seg3.sel.sp@data$Class != 1 & seg3.sel.sp@data$Class != 30 & seg3.sel.sp@data$Class != 2 & seg3.sel.sp@data$Class != 7 & seg3.sel.sp@data$Cl1Fl == 1)] <- 333

# Get border of flow objects to class 1 and 30
L3.BodyFlow.RelClass <- RelationalClassFunction(spdf = seg3.sel.sp.body, nb =  L3.nb, class.var = 130, class.var2 = 333, var = c("DTM", "Area"), quiet = FALSE, class.as.neighbors = FALSE)
L3.BodyFlow.RelClass.1.30 <- L3.BodyFlow.RelClass[c("ID", "bor_cl_rel")]
colnames(L3.BodyFlow.RelClass.1.30) <- c("ID", "rB_Fl_130")


# merge to table
seg3.sel.sp@data <- merge(seg3.sel.sp@data, L3.BodyFlow.RelClass.1.30 , by = "ID", all.x = TRUE)




# ... Set Class 3! --------------------------------------------------
# "Cl1Fl" = 1 AND "rb_Fl_130" <> -9999 AND ("rb_Fl_130"  >= 60 
# OR "rBor_Cl10"  <> "rB_Fl_130" AND "Cl1FlInv" = 0 AND   "Sh_Ind" < 4.472944 AND "rB_Fl_130" >=  51 AND "Cl1FlNb" = 1 AND "Slp" >= 8.821035)

# "rB_Fl_130" >= 20 AND "SfRghn" >= 0.01109138 AND "Sh_Ind" <= 3.453671 AND "Slp"  >= 8.821035 AND "Slp"  < 23.176670 AND "Cl1FlInv" = 0
# OR "rB_Fl_130" >= 20 AND "Slp_SD" >= 4.919841 AND "Slp"  >= 8.821035 AND "Slp"  < 23.176670 AND "Cl1FlInv" = 0
# OR "rB_Fl_130" >= 20 AND "Area" > 4300 AND "Sh_Ind" <= 3.453671 AND "Cl1FlNb" = 1

seg3.sel.sp$Class[which(seg3.sel.sp@data$Class != 1 & seg3.sel.sp@data$Class != 30 & seg3.sel.sp@data$Class != 2 & seg3.sel.sp@data$Class != 7 & seg3.sel.sp@data$rB_Fl_130 > 0 
                        & ((seg3.sel.sp@data$Cl1Fl == 1 & seg3.sel.sp@data$rB_Fl_130 >= 60)
                           | (seg3.sel.sp@data$rBor_Cl10 != seg3.sel.sp@data$rB_Fl_130 & seg3.sel.sp@data$Cl1FlInv == 0 & seg3.sel.sp@data$Sh_Ind < 4.472944 & seg3.sel.sp@data$rB_Fl_130 >= 51 & seg3.sel.sp@data$Cl1FlNb == 1 & seg3.sel.sp@data$Slp >= 8.821035)
                           
                           | (seg3.sel.sp@data$rB_Fl_130 >= 20 &  seg3.sel.sp@data$Slp >= 8.821035 & seg3.sel.sp@data$Slp < 23.176670 & ((seg3.sel.sp@data$SfRghn >= 0.01109138 & seg3.sel.sp@data$Sh_Ind <= 3.453671 & seg3.sel.sp@data$Cl1FlInv == 0) 
                                                                                                                                         | (seg3.sel.sp@data$Slp_SD >= 4.919841 & seg3.sel.sp@data$Cl1FlInv == 0)))
                           | (seg3.sel.sp@data$rB_Fl_130 >= 20 & seg3.sel.sp@data$Area >= 4300 & seg3.sel.sp@data$Sh_Ind <= 3.453671 & seg3.sel.sp@data$Cl1FlNb == 1)
                        ))] <- 3



seg3.sel.sp$Class[which(seg3.sel.sp@data$Class == 30)] <- 3







# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# L3: Neighbor Growing and final Cleaning ------------------------------------

# ... Neighbor Growing --------------------------------------------
# growing of flanks (class 22)
L3.NeighborGrowing.LS <- NeighborGrowing(spdf = subset(seg3.sel.sp, seg3.sel.sp@data$Class == 1 | seg3.sel.sp@data$Class == 2 | seg3.sel.sp@data$Class == 3), 
                                         ID.start = seg3.sel.sp@data$ID[which(seg3.sel.sp@data$Class == 1)], return.gUnaryUnionNeighbors = TRUE, return.input = FALSE)


# necessairy to write growing data!
rgdal::writeOGR(L3.NeighborGrowing.LS, dsn = paste(getwd(), dirname(L3.landslide), sep = "/"), overwrite_layer = TRUE, 
                layer = file_path_sans_ext(basename(L3.landslide)), driver = "ESRI Shapefile")




# shape indices
rsaga.geoprocessor(lib="shapes_polygons", module = 7, env = env, show.output.on.console = FALSE, param = list(
  SHAPES = L3.landslide, INDEX = L3.landslide))



correctDBF(L3.landslide, adjust.n = 8, new.colnames = c("Area", "P","P_A", "P_sqrt_A", "Mx_Dist","D_A", "D_sqrt_A", "Sh_Ind"))









# OUT -----------------
# # #
# # transform NA to -9999
seg3.sel.sp.out <- seg3.sel.sp
replaceInvalids(seg3.sel.sp.out@data, replace.value = -9999)

rgdal::writeOGR(seg3.sel.sp.out, dsn = paste(getwd(), dirname(seg3.sel), sep = "/"), overwrite_layer = TRUE,
                layer = "L3_seg_sel_out", driver = "ESRI Shapefile")
# correctDBF(paste0(file_path_sans_ext(seg3.sel), "_out.shp"), end.n = 0, new.colnames = colnames(seg3.sel.sp@data))
# # #








# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 3 Accuracy Assessment ------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

AccAss.flank <- read.dbf("Output/Accuracy_Assessment/AccAss_flank.dbf")$dbf
AccAss.scarp <- read.dbf("Output/Accuracy_Assessment/AccAss_scarp.dbf")$dbf
AccAss.landslide <- read.dbf("Output/Accuracy_Assessment/AccAss_landslide.dbf")$dbf
AccAss.ratio <- read.dbf("Output/Accuracy_Assessment/AccAss_ratio.dbf")$dbf


AccAss.inventory <- read.dbf("Output/Accuracy_Assessment/inventory.dbf")$dbf


# # # comparison true false of detection
AccAss.landslide.true <- data.frame(Area = AccAss.landslide$Area[which(AccAss.landslide$AccAss == 1)])
AccAss.landslide.true$AccAss <- "true"

AccAss.landslide.false <- data.frame(Area = AccAss.landslide$Area[which(AccAss.landslide$AccAss == 0)])
AccAss.landslide.false$AccAss <- "false"


AccAss.landslide.true_false <- rbind(AccAss.landslide.true, AccAss.landslide.false)

# http://stackoverflow.com/questions/3541713/how-to-plot-two-histograms-together-in-r
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# http://docs.ggplot2.org/current/guide_legend.html
ggplot(AccAss.landslide.true_false, aes(Area, fill = AccAss)) + geom_density(alpha = 0.6) +
      theme(axis.title=element_text(face="bold", size="25"), axis.text=element_text(size=23), legend.position="bottom", legend.box="horizontal", legend.title = element_text(size = 30, face="bold"), legend.text=element_text(size=26)) + 
            scale_fill_manual(values=c("#FF3333", "#00CC00"), name="detected landslide", labels = c("false  ", "true")) + xlab("Area (m²)") +
            guides(fill = guide_legend(title.position = "top"))



# # # comparison of true and inventory
AccAss.inventory <- data.frame(Area = AccAss.inventory$Area)
AccAss.inventory$AccAss <- "inventory"


AccAss.landslide.true_inventory <- rbind(AccAss.inventory, AccAss.landslide.true)

# # # ggplots2

ggplot(AccAss.landslide.true_inventory, aes(Area, fill = AccAss)) + geom_density(alpha = 0.4, adjust = 4) +
  theme(axis.title=element_text(face="bold", size="25"), axis.text=element_text(size=23), legend.direction = "horizontal", legend.position="bottom", legend.title = element_text(size = 30, face="bold"), legend.text=element_text(size=26)) + 
  scale_fill_manual(values=c("#CC9900", "#00CC00"), name="", labels = c("landslides of inventory  ", "correctly classified landslides")) + xlab("Area (m²)")

# true: correctly classified landslides ... remove title

# # # comparison of true, false and inventory
AccAss.landslide.true_false_inventory <- rbind(AccAss.inventory, AccAss.landslide.true, AccAss.landslide.false)

ggplot(AccAss.landslide.true_false_inventory, aes(Area, fill = AccAss)) + geom_density(alpha = 0.4, adjust = 1) +
  theme(axis.title=element_text(face="bold", size="11"), legend.position="right", legend.title = element_text(size = 11, face="bold")) + 
  scale_fill_manual(values=c("#CC9900", "#00CC00", "#FF3333"), name="Detected\nlandslide") + xlab("Area (m²)")





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Save Data ---------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
print("... Save Data")
save.image("L3.RData")
