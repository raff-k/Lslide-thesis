# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# THESIS: FIRST SEGMENTATION LEVEL - PART 1 ------------------------------------------------------------------------
#
# Raphael Knevels
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# DESCRIPTION:
# This is the initial script for the automated, object-oriented detection
# of landslides. Here, the following parameters are computed:
#   - initialization of variables and software (i.e. GRASS GIS settings)
#   - land surface parameters, textural features, and shape metrics
#   - fine-scale segmentation in L1
#   - zonal statistic for objects resulted from segmentation L1
#   - k-means thresholds for the further analysis
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Date: 20.10.2016 - 20.04.2017
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# CONTENT -----------------------------------------------------------------
# 1 PACKAGES, FUNCTIONS & VARIABLES
# 2 CALCULATION OF PARAMETERS
# 3 L1: OBJECT-ORIENTED IMAGE ANALYSIS
#   ... grid statistics
#   ... selection based on maximum curvature
#   ... compute main direction and flow direction
#   ... k-means thresholding






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 1 PACKAGES, FUNCTIONS & VARIABLES ------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, maptools, mclust, raster, rgdal, rgeos, rgrass7, RSAGA, tools, sf, sp, shapefiles)

# path
setwd("/home/raphael/Data")



# ... variables --------------------------------------------------------
print("Initialize Variables")

dtm <- "Output/DTM/dtm.sgrd"
dtm.GRASS <- paste0(file_path_sans_ext(dtm), ".sdat")
dtm.tif <- "Output/DTM/dtm.tif"

dtm51 <- "Output/DTM/dtm_51.sgrd"
dtm.dif.dtm51 <- "Output/DTM/dtm_dif_dtm51.sgrd"
dtm.D51.edge.vigra <- "Output/DTM/dtm_D51_edge_vigra.sgrd"
dtm.D51.edge.vigra.re03 <- "Output/DTM/dtm_D51_edge_vigra_re03.sgrd"
dtm.texture <- "Output/DTM/dtm_texture.sgrd"

slope <- "Output/Slope/slope.sgrd"
slope.GRASS <- paste0(file_path_sans_ext(slope) , ".sdat")
slope.edge.vigra <- "Output/Slope/slope_edge_vigra.sgrd"
slope.edge.vigra.re03 <- "Output/Slope/slope_edge_vigra_re03.sgrd"
slope.texture <- "Output/Slope/slope_texture.sgrd"
slp.tif <- "Output/Slope/slope.tif"
texture.slp.num.quantile <- round(maxValue(raster(slp.tif)) - minValue(raster(slp.tif)))


curv.max15 <- "Output/Other/curv_max15.sgrd"
curv.max15.IQR<- "Output/Other/curv_max15_IQR.sgrd"
curv.max15.IQR.buf <- "Output/Other/curv_max15_IQR_buf.sgrd"
curv.max15.IQR.buf.dist <- "5"
curv.max15.IQR.buf.shp <- "Output/Other/curv_max15_IQR_buf.shp"

# need TauDEM 5.x, no space inside variable names!
dtm.PitRemove <- "Output/DTM/dtm_fill.tif"
flow.Dinf.rad <- "Output/Other/flow_Dinf_rad.tif"
flow.Dinf.deg <- "Output/Other/flow_Dinf_deg.tif"
flow.Dinf.Area <- "Output/Other/flow_Dinf_Area.tif"
flow.Dinf.sin <- "Output/Other/flow_Dinf_sin.tif"
flow.Dinf.cos <- "Output/Other/flow_Dinf_cos.tif"


surf.roughn15 <- "Output/Other/surf_roughn15.sgrd"
surf.roughn3 <- "Output/Other/surf_roughn3.sgrd"


segmentation.first <- "Output/Segmentation/L1_seg.shp"
segmentation.first.grid <- "Output/Segmentation/L1_seg.sgrd"
seg1.sel <- "Output/Segmentation/L1_seg_sel.shp"

segmentation.second <- "Output/Segmentation/L2_seg.shp"
segmentation.second.grid <- "Output/Segmentation/L2_seg.sgrd"
seg2.sel <- "Output/Segmentation/L2_seg_sel.shp"

texture.entr.flow <- paste0(getwd(), "/", "Output/Other")
texture.entr.win <- 3
texture.entr.dist <- 1
texture.sa.flow <- paste0(getwd(), "/", "Output/Other")
texture.sa.win <- 5
texture.sa.dist <- 1
texture.output.name <- c("textureFlowDir", "textureFlowDirPer") # same as default
# texture.dtm.num.quantile <- round(maxValue(raster(dtm.tif)) - minValue(raster(dtm.tif)))
texture.dtm.ruleset <- "Output/DTM/dtm_texture_ruleset.txt"
texture.slp.ruleset <- "Output/Slope/slp_texture_ruleset.txt"
texture.calc <- TRUE
texture.entr <- paste0(texture.entr.flow, "/",  texture.output.name[1], texture.entr.win, texture.sa.dist, "_Entr.tiff")
texture.sa <- paste0(texture.sa.flow, "/",  texture.output.name[1], texture.sa.win, texture.sa.dist, "_SA.tiff")


skyviewfactor <- paste0(getwd(), "/", "Output/DTM/dtm_SVF_R20_D16.tif")



L1.final.grid <- "Output/Segmentation/L1_final.sgrd"
L1.final.grid.buf <- "Output/Segmentation/L1_final_buf.sgrd"
L1.final.grid.buf.dist <- 5
L1.final <- "Output/Segmentation/L1_final.shp" 




# .... initialize SAGA & GRASS --------------------------------------------------------
print("Initialize SAGA and GRASS")

# SAGA
env <-   rsaga.env()
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
    if(nrow(m) > size.sample)
    {
      set.seed(seed)
      #        User      System verstrichen 
      #       277.97        1.23      281.60 
      n.clust <- mclust::Mclust(m, G = 1:G, initialization = list(subset = sample(1:nrow(m), size = size.sample)))
    } else
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
  # found here: http://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table
  
  # NA and NAN
  for(j in seq_len(ncol(x)))
    set(x, which(is.na(x[[j]])), j, replace.value)
  
  # NULL
  for(j in seq_len(ncol(x)))
    set(x, which(is.null(x[[j]])), j, replace.value)
  
}









# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 2 CALCULATION OF PARAMETERS ------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# ... openness as difference between dtm and dtm50 by Van den Eeckhaut et al. (2012) ---------------------------------------------
### dtm 50 - CAN TAKE A WHILE!
# rsaga.get.usage("grid_filter", 0, env = env)
# METHOD: [0] Smooth
# "------ Run of Segmentation: 27.934 Minutes ------"
if(!file.exists(file.path(getwd(), dtm51)))
{
  # print(parseGRASS("r.neighbors"))
  execGRASS('r.neighbors', flags = c("quiet", "overwrite"), parameters = list(
    input = "dtm", size = 51, output = "dtm51", method = "average"))
  
  dtm51.GRASS <- paste0(getwd(), "/", file_path_sans_ext(dtm51), ".sdat")
  
  # print(parseGRASS("r.out.gdal"))
  execGRASS('r.out.gdal',  flags=c("overwrite", "quiet"), parameters = list(
    input = "dtm51", output = dtm51.GRASS, format = "SAGA"))
  
  rm(dtm51.GRASS)
}

if(!file.exists(file.path(getwd(), dtm.dif.dtm51)))
{
  # formula for calculation
  formula.expression <- paste0("a-b")
  
  # rsaga.get.usage("grid_calculus", 1, env = env)
  rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = FALSE, param=list(
    GRIDS = paste(c(dtm, dtm51), collapse = ";"), RESULT = dtm.dif.dtm51, FORMULA = formula.expression, FNAME = "1"))
  
  rm(formula.expression)
}


## edge detection (vigra)
# rsaga.get.modules("imagery_vigra", env = env)
# rsaga.get.usage("imagery_vigra", 1, env = env)
# TYPE: [0] Canny | [1] Shen-Castan
if(!file.exists(file.path(getwd(), dtm.D51.edge.vigra)))
{
  rsaga.geoprocessor(lib="imagery_vigra", module = 1, env = env, show.output.on.console = FALSE, param = list(
    INPUT = dtm.dif.dtm51, OUTPUT = dtm.D51.edge.vigra, TYPE = "0", SCALE = "1.0", THRESHOLD = "0.05"))
}


## reclassify vigra output for GRASS, because GRASS cannot handle No-Data in segmentation
# rsaga.get.modules("grid_tools", env = env)
# rsaga.get.usage("grid_tools", 15, env = env)
# METHOD: [0] single | SOPERATOR: [0] =
if(!file.exists(file.path(getwd(), dtm.D51.edge.vigra.re03)))
{
  rsaga.geoprocessor(lib="grid_tools", module = 15, env = env, show.output.on.console = FALSE, param = list(
    INPUT = dtm.D51.edge.vigra, RESULT = dtm.D51.edge.vigra.re03, METHOD = "0", OLD = "1.0", NEW = "3.0",
    SOPERATOR = "0", NODATAOPT = "1", NODATA = "0"))
}






# ... slope and slope derivatives ---------------------------------------------
## slope
# rsaga.get.usage("ta_morphometry", 0, env = env)
# Method 6: [6] 9 parameter 2nd order polynom (Zevenbergen & Thorne 1987) 
# Unit is [1] degree
if(!file.exists(file.path(getwd(), slope)))
{
  rsaga.geoprocessor(lib="ta_morphometry", env = env, module=0, show.output.on.console = FALSE, param=list(
    ELEVATION= dtm, SLOPE = slope, METHOD = "6", UNIT_SLOPE = "1"))
}

## edge detection (vigra)
# rsaga.get.modules("imagery_vigra", env = env)
# rsaga.get.usage("imagery_vigra", 1, env = env)
# TYPE: [0] Canny | [1] Shen-Castan
if(!file.exists(file.path(getwd(), slope.edge.vigra)))
{
  rsaga.geoprocessor(lib="imagery_vigra", module = 1, env = env, show.output.on.console = FALSE, param = list(
    INPUT = slope, OUTPUT = slope.edge.vigra, TYPE = "0", SCALE = "1.0", THRESHOLD = "1.0"))
}

## reclassify vigra output for GRASS, because GRASS cannot handle No-Data in segmentation
# rsaga.get.modules("grid_tools", env = env)
# rsaga.get.usage("grid_tools", 15, env = env)
# METHOD: [0] single | SOPERATOR: [0] =

if(!file.exists(file.path(getwd(), slope.edge.vigra.re03)))
{
  rsaga.geoprocessor(lib="grid_tools", module = 15, env = env, show.output.on.console = FALSE, param = list(
    INPUT = slope.edge.vigra, RESULT = slope.edge.vigra.re03, METHOD = "0", OLD = "1.0", NEW = "3.0",
    SOPERATOR = "0", NODATAOPT = "1", NODATA = "0"))
}



# ... D-infinity flow ---------------------------------------------
## fill DTM
if(!file.exists(file.path(getwd(), dtm.PitRemove)))
{
  command <- paste("mpiexec PitRemove -z", paste(getwd(), dtm.tif, sep = "/"), "-fel", paste(getwd(), dtm.PitRemove, sep = "/"))
  system(command, show.output.on.console = FALSE)
}

## calculate D-Infinity Flow Direction
if(!file.exists(file.path(getwd(), flow.Dinf.rad)))
{
  command <- paste("mpiexec DinfFlowDir -fel", paste(getwd(), dtm.PitRemove, sep = "/"), "-ang", paste(getwd(), flow.Dinf.rad, sep = "/"))
  system(command, show.output.on.console = FALSE)
  
  # read raster
  # flow.Dinf <- raster::raster(paste(getwd(), flow.Dinf.rad, sep = "/"))
  # proj4string(flow.Dinf) <- CRS(paste0("+init=", epsg.code))
  flow.Dinf <- rgdal::readGDAL(paste(getwd(), flow.Dinf.rad, sep = "/"), silent = TRUE)
  
  
  # create cosine and sinus raster for mean value
  # sin.flow.Dinf <- sin(flow.Dinf)
  # cos.flow.Dinf <- cos(flow.Dinf)
  sin.flow.Dinf <- flow.Dinf
  sin.flow.Dinf@data$band1 <- sin(sin.flow.Dinf@data$band1)
  cos.flow.Dinf <- flow.Dinf
  cos.flow.Dinf@data$band1 <- cos(cos.flow.Dinf@data$band1)
  
  
  # writeRaster(sin.flow.Dinf, flow.Dinf.sin, overwrite = TRUE)
  # writeRaster(cos.flow.Dinf, flow.Dinf.cos, overwrite = TRUE)
  rgdal::writeGDAL(dataset = sin.flow.Dinf, fname = paste0(getwd(), "/", flow.Dinf.sin))
  rgdal::writeGDAL(dataset = cos.flow.Dinf, fname = paste0(getwd(), "/", flow.Dinf.cos))
  
  # remove data from environment
  rm(flow.Dinf, sin.flow.Dinf, cos.flow.Dinf )
}

# transformation from radiant to degree AND switch direction
if(!file.exists(file.path(getwd(), flow.Dinf.deg)))
{
  # read raster
  # flow.Dinf <- raster::raster(paste(getwd(), flow.Dinf.rad, sep = "/"))
  flow.Dinf <- rgdal::readGDAL(paste(getwd(), flow.Dinf.rad, sep = "/"), silent = TRUE)
  # proj4string(flow.Dinf) <- CRS(paste0("+init=", epsg.code))
  
  # rotate and write raster
  flow.Dinf@data$band1 <- ((flow.Dinf@data$band1 * (-180)/pi) + 90) %% 360 # from counter-clockwise to clockwise AND from E = 0 to N = 0 degree
  rgdal::writeGDAL(dataset = flow.Dinf, fname = paste0(getwd(), "/", flow.Dinf.deg))
  # rgdal::writeGDAL(dataset = flow.Dinf, fname = flow.Dinf.deg)
  
  # writeRaster(flow.Dinf.Rotated, flow.Dinf.deg, overwrite = TRUE)
  
  # remove data from environment
  rm(flow.Dinf)
}


# calculate watershed | contributing area
if(!file.exists(file.path(getwd(), flow.Dinf.Area)))
{
  command <- paste("mpiexec AreaDinf -ang", paste(getwd(), flow.Dinf.rad, sep = "/"), "-sca", paste(getwd(), flow.Dinf.Area, sep = "/"))
  system(command, show.output.on.console = FALSE)
}




# ... Curvature ---------------------------------------------
## Curvature Max by Wood (1996) - paramteres by Tarolli (2012)
if(!file.exists(file.path(getwd(), curv.max15)))
{
  # rsaga.get.usage("ta_morphometry", 23, env = env)
  rsaga.geoprocessor(lib = "ta_morphometry", env = env, module = 23, show.output.on.console = FALSE, param=list(
    DEM = dtm, MAXIC = curv.max15, SIZE = 15))
}

if(!exists("curvature.max15.IQR.threshold")) # 0.004074238
{
  # read max curvature file
  curv.max15.r <- rgdal::readGDAL(paste0(file_path_sans_ext(curv.max15), ".sdat"), silent = TRUE)
  curv.max15.stat<- quantile(data.table(curv.max15.r@data$band1), na.rm = TRUE)
  
  # calculate threshord by interquantile range
  curv.max15.IQR.m <- 1.5
  curv.max15.IQR.thr <- curv.max15.IQR.m *(curv.max15.stat[[4]] - curv.max15.stat[[2]]) # 3. third - first quartile # 0.007163971
  
  rm(curv.max15.r, curv.max15.stat, curv.max15.IQR.m) # clean up environment
}

# calculate true/false for maximal curvature (15x15) by IQR.threshold
if(!file.exists(file.path(getwd(), curv.max15.IQR)))
{
  # formula for calculation
  formula.expression <- paste0("ifelse(gt(a,", curv.max15.IQR.thr, "), 1, (-99999))")
  #formula.expression <- paste0("gt(g1,",  curvature.max15.IQR.threshold, ")")
  
  # rsaga.get.usage("grid_calculus", 1, env = env)
  rsaga.geoprocessor(lib = "grid_calculus", env = env, module = 1, show.output.on.console = FALSE, param=list(
    GRIDS = curv.max15, RESULT = curv.max15.IQR, FORMULA = formula.expression, FNAME = "1"))
  
  rm(formula.expression)
}


# convert buffered maximal curvature (15x15) by IQR.threshold to shapefile
if(!file.exists(file.path(getwd(), curv.max15.IQR.buf.shp)))
{
  # grid buffering
  # rsaga.get.usage("grid_tools", 8, env = env)
  # BUFFERTYPE: [0] Fixed | [1] Cell value 
  rsaga.geoprocessor(lib="grid_tools", module = 8, env = env, show.output.on.console = FALSE, param = list(
    FEATURES = curv.max15.IQR, BUFFER = curv.max15.IQR.buf, DIST = curv.max15.IQR.buf.dist, BUFFERTYPE = "0")) 
  
  
  # rsaga.get.usage("shapes_grid", 6, env = env)
  rsaga.geoprocessor(lib="shapes_grid", module = 6, env = env, show.output.on.console = FALSE, param = list(
    GRID = curv.max15.IQR.buf, POLYGONS = curv.max15.IQR.buf.shp)) 
  
  # correct corrupt field naming of SAGA GIS
  correctDBF(x = curv.max15.IQR.buf.shp, end.n = 0, new.colnames = c("Curv_M", "ID", "NAME"))
  
}



# ... Surface Roughness ---------------------------------------------
# Grass GIS Add-on: r.vector.ruggedness - Vector Ruggedness Measure by Sappington et al. (2007) 
# print(parseGRASS("g.extension"))
# execGRASS('g.extension', flags = c("quiet"), extension = 'r.vector.ruggedness')

if(!file.exists(file.path(getwd(), surf.roughn15)))
{
  # print(parseGRASS("r.vector.ruggedness"))
  execGRASS('r.vector.ruggedness', flags = c("quiet", "overwrite"), parameters = list(
    elevation = "dtm", size = 15, output = "ruggedness15"))
  
  surf.roughn15.GRASS <- paste0(getwd(), "/", file_path_sans_ext(surf.roughn15), ".sdat")
  
  # print(parseGRASS("r.out.gdal"))
  execGRASS('r.out.gdal',  flags=c("overwrite", "quiet"), parameters = list(
    input = "ruggedness15", output = surf.roughn15.GRASS, format = "SAGA"))
  
  rm(surf.roughn15.GRASS)
}


if(!file.exists(file.path(getwd(), surf.roughn3)))
{
  # print(parseGRASS("r.vector.ruggedness"))
  execGRASS('r.vector.ruggedness', flags = c("quiet", "overwrite"), parameters = list(
    elevation = "dtm", size = 3, output = "ruggedness7"))
  
  surf.roughn3.GRASS <- paste0(getwd(), "/", file_path_sans_ext(surf.roughn3), ".sdat")
  
  # print(parseGRASS("r.out.gdal"))
  execGRASS('r.out.gdal',  flags=c("overwrite", "quiet"), parameters = list(
    input = "ruggedness7", output = surf.roughn3.GRASS, format = "SAGA"))
  
  rm(surf.roughn3.GRASS)
}



# ... SkyViewFactor ---------------------------------------------
skyViewFactor(path.software = "E:/Masterarbeit/Software/SkyViewFactor",  sky_view_factor = c(1, 16, 20, 0, "low"),
              path.input = "E:/Masterarbeit/Data/Output/DTM/dtm.tif", quiet = FALSE)



# ... Texture ---------------------------------------------
# the result of the texture is affected on the min-max value of the input dtm.
# that means that textures of a representative region differ from that of a whole extent.
# to avoid this, the whole dtm should be first "recoded" by a ruleset, which is then also used to recode the region dtm

# create ruleset (should be done for whole extent)
if(!file.exists(file.path(getwd(), texture.dtm.ruleset)))
{
  # print(parseGRASS("r.quantile"))
  texture.dtm.quantile <- execGRASS('r.quantile', flags = c("r", "quiet"), intern = TRUE, parameters = list(
    input = "dtm", quantiles = texture.dtm.num.quantile))
  
  
  write(texture.dtm.quantile, texture.dtm.ruleset)
  
  rm(texture.dtm.quantile)
}  

if(!file.exists(file.path(getwd(), dtm.texture)))
{    
  texture.dtm.ruleset.GRASS <- paste0(getwd(), "/", texture.dtm.ruleset)
  dtm.texture.GRASS <- paste0(getwd(), "/", file_path_sans_ext(dtm.texture), ".sdat")
  
  # print(parseGRASS("r.recode"))
  execGRASS('r.recode', flags = c("overwrite", "d"), parameters = list(
    input = "dtm", output = "dtm_texture", rules = texture.dtm.ruleset.GRASS))
  
  
  execGRASS('r.out.gdal',  flags=c("overwrite", "quiet"), parameters = list(
    input = "dtm_texture", output = dtm.texture.GRASS, format = "SAGA"))
  
  rm(texture.dtm.ruleset.GRASS, dtm.texture.GRASS)
}



if(texture.calc == TRUE)
{
  # Entropy
  TextureFlow(grass.texture.flowDir = paste0(getwd(), "/", flow.Dinf.deg), grass.texture.input = "dtm_texture", grass.texture.method = "entr", grass.texture.name = c("_Entr_"),
              grass.texture.window = texture.entr.win, grass.texture.distance = texture.entr.dist, grass.texture.save = TRUE, grass.texture.save.path = texture.entr.flow, texture.output.name = texture.output.name)
  
  
  
  # Sum Average
  TextureFlow(grass.texture.flowDir = paste0(getwd(), "/", flow.Dinf.deg), grass.texture.input = "dtm_texture", grass.texture.method = "sa", grass.texture.name = c("_SA_"),
              grass.texture.window = texture.sa.win, grass.texture.distance = texture.sa.dist, grass.texture.save = TRUE, grass.texture.save.path = texture.sa.flow, texture.output.name = texture.output.name)
}  






# ... Texture SLOPE ---------------------------------------------
# the result of the texture is affected on the min-max value of the input dtm.
# that means that textures of a representative region differ from that of a whole extent.
# to avoid this, the whole dtm should be first "recoded" by a ruleset, which is then also used to recode the region dtm

# create ruleset (should be done for whole extent)
if(!file.exists(file.path(getwd(), texture.slp.ruleset)))
{
  # print(parseGRASS("r.quantile"))
  texture.slp.quantile <- execGRASS('r.quantile', flags = c("r", "quiet"), intern = TRUE, parameters = list(
    input = "slope", quantiles = texture.slp.num.quantile))
  
  
  write(texture.slp.quantile, texture.slp.ruleset)
  
  rm(texture.slp.quantile)
}  

if(!file.exists(file.path(getwd(), slp.texture)))
{    
  texture.slp.ruleset.GRASS <- paste0(getwd(), "/", texture.slp.ruleset)
  slp.texture.GRASS <- paste0(getwd(), "/", file_path_sans_ext(slope.texture), ".sdat")
  
  # print(parseGRASS("r.recode"))
  execGRASS('r.recode', flags = c("overwrite", "d"), parameters = list(
    input = "slope", output = "slope_texture", rules = texture.slp.ruleset.GRASS))
  
  
  execGRASS('r.out.gdal',  flags=c("overwrite", "quiet"), parameters = list(
    input = "slope_texture", output = slp.texture.GRASS, format = "SAGA"))
  
  rm(texture.slp.ruleset.GRASS, slp.texture.GRASS)
}



if(texture.calc == TRUE)
{
  # Entropy
  TextureFlow(grass.texture.flowDir = paste0(getwd(), "/", flow.Dinf.deg), grass.texture.input = "dtm_texture", grass.texture.method = "entr", grass.texture.name = c("_Entr_"),
              grass.texture.window = texture.entr.win, grass.texture.distance = texture.entr.dist, grass.texture.save = TRUE, grass.texture.save.path = texture.entr.flow, texture.output.name = texture.output.name)
  
  
  
  # Sum Average
  TextureFlow(grass.texture.flowDir = paste0(getwd(), "/", flow.Dinf.deg), grass.texture.input = "dtm_texture", grass.texture.method = "sa", grass.texture.name = c("_SA_"),
              grass.texture.window = texture.sa.win, grass.texture.distance = texture.sa.dist, grass.texture.save = TRUE, grass.texture.save.path = texture.sa.flow, texture.output.name = texture.output.name)
}  






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 3 OBJECT-ORIENTED IMAGE ANALYSIS ------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ... SEGMENTATION LEVEL 1 ------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

print("Start SEGMENTATION LEVEL 1 ...")

# ... L1: perform segmentation for level 1 ----------------------------------------
print("... Segmentation")
if(!file.exists(file.path(getwd(), segmentation.first)))
{
  segmentation(Tool = "SAGA", Segments.Grid = segmentation.first.grid, Segments.Poly = segmentation.first, Input.Grid = c(slope, slope.edge.vigra),
               Seed.Method = "Fast Representativeness", Fast.Representativeness.LevelOfGeneralisation = 1.25,  Saga.Segmentation.Method = "0", Saga.Segmentation.Sig.2 = "125",
               Generalisation.Flac = FALSE, Saga.Segmentation.Leafsize = 1024, NoData = TRUE, Mask = slope)
}






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ... L1: calculate grid statistics for curvature max IQR ----------------------------------------
print("... Curvature Max IQR")
#rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
#  GRIDS = curv.max15.IQR.buf, POLYGONS = segmentation.first, METHOD = "0", NAMING = "1", RESULT = segmentation.first,
#  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "1", MEAN = "0", VAR = "0", STDDEV = "0", QUANTILE = 0))


#correctDBF(segmentation.first, adjust.n = 1, new.colnames = c("Curv_Sum"))






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ... L1: pre-selection of first segmentation for speed-up and memory-reasons ----------------------------------------
# load result of first segmentation
print("... Selection of Polygons that intersect Curvature Max IQR")
#if(!exists("seg1.sf"))
#{
  seg1.sf <- sf::st_read(dsn = paste(getwd(), dirname(segmentation.first), sep = "/"), layer = file_path_sans_ext(basename(segmentation.first)), quiet = TRUE, stringsAsFactors = FALSE)
  seg1.sf <- sf::st_transform(seg1.sf, sp::CRS(paste0("+init=", tolower(epsg.code)))@projargs)
  
   # transform to SP
   seg1.sel.sp <- as(seg1.sf, 'Spatial')
#}

#if(!file.exists(file.path(getwd(), seg1.sel)))
#{
   print("... selection")
  # select all polygon that intersected with the max curvature raster
  seg1.sel.sp <-  subset(seg1.sel.sp, seg1.sel.sp@data$Curv_Sum > 0)

  rgdal::writeOGR(seg1.sel.sp, dsn = paste(getwd(), dirname(seg1.sel), sep = "/"), overwrite_layer = TRUE, 
                  layer = file_path_sans_ext(basename(seg1.sel)), driver = "ESRI Shapefile")
  
  print("... remove objects")
  rm(seg1.sf,seg1.sel.sp) # ,segmnetation.frist.sf.selection)
#}







# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ... L1: calculate grid statistics for slope, dtm, and flow direction ----------------------------------------
print("... Calculation of Grid Statistics")


# slope
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = slope, POLYGONS = seg1.sel, METHOD = "0", NAMING = "1", RESULT = seg1.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))



# CurvMax 
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = curv.max15, POLYGONS = seg1.sel, METHOD = "0", NAMING = "1", RESULT = seg1.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# dtm
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = dtm, POLYGONS = seg1.sel, METHOD = "0", NAMING = "1", RESULT = seg1.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# dtm dif 51
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = dtm.dif.dtm51, POLYGONS = seg1.sel, METHOD = "0", NAMING = "1", RESULT = seg1.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# flow direction
# flow.Dinf.cos <- paste0(file_path_sans_ext(flow.Dinf.cos), ".sgrd")
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = flow.Dinf.cos, POLYGONS = seg1.sel, METHOD = "0", NAMING = "1", RESULT = seg1.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))

# flow.Dinf.sin <- paste0(file_path_sans_ext(flow.Dinf.sin), ".sgrd")
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = flow.Dinf.sin, POLYGONS = seg1.sel, METHOD = "0", NAMING = "1", RESULT = seg1.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# vector ruggedness measure 
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = surf.roughn15, POLYGONS = seg1.sel, METHOD = "0", NAMING = "1", RESULT = seg1.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# texture entropy
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = texture.entr, POLYGONS = seg1.sel, METHOD = "0", NAMING = "1", RESULT = seg1.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# texture sa
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = texture.sa, POLYGONS = seg1.sel, METHOD = "0", NAMING = "1", RESULT = seg1.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# sky view factor
rsaga.geoprocessor(lib="shapes_grid", module = 2, env = env, show.output.on.console = FALSE, param = list(
  GRIDS = skyviewfactor, POLYGONS = seg1.sel, METHOD = "0", NAMING = "1", RESULT = seg1.sel,
  COUNT = "0", MIN = "0", MAX = "0", RANGE = "0", SUM = "0", MEAN = "1", VAR = "0", STDDEV = "0", QUANTILE = 0))


# shape indices
rsaga.geoprocessor(lib="shapes_polygons", module = 7, env = env, show.output.on.console = FALSE, param = list(
  SHAPES = seg1.sel, INDEX = seg1.sel))


correctDBF(seg1.sel, adjust.n = 18, new.colnames = c("Slp", "CurvM15","DTM", "DTM_D_51", "Fl_Cos", "Fl_Sin", "SfRghn", "t_Entr", "t_SA", "SVF",
                                                     "Area", "P","P_A", "P_sqrt_A", "Mx_Dist", "D_A", "D_sqrt_A", "Sh_Ind"))






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ... L1: main direction and flow direction ----------------------------------------
print("... Computation of main flow and mean direction")

# re-load polygons
# rm(seg1.sel.sf)
seg1.sel.sf <- sf::st_read(dsn = paste(getwd(), dirname(seg1.sel), sep = "/"),
                           layer = file_path_sans_ext(basename(seg1.sel)), quiet = TRUE)

# convert simple feature to spatial
seg1.sel.sp <- as(seg1.sel.sf, "Spatial")

# set NAs
seg1.sel.sp@data$Fl_Cos[((seg1.sel.sp@data$Fl_Cos == 0) & (seg1.sel.sp@data$Fl_Sin == 0)) | ((seg1.sel.sp@data$Fl_Cos == -9999) & (seg1.sel.sp@data$Fl_Sin == -9999))] <- NA
seg1.sel.sp@data$Fl_Sin[is.na(seg1.sel.sp@data$Fl_Cos)] <- NA

# calculate mean Flow
seg1.sel.sp@data$Flow <- ((atan2(seg1.sel.sp@data$Fl_Sin, seg1.sel.sp@data$Fl_Cos) * (-180)/pi) + 90) %% 360



# caluclation of object orientation and objects lying perpendicular to flow direction
seg1.sel.MainDir <- MainDirection(seg1.sel.sp)
seg1.sel.sp@data$MnDir <- seg1.sel.MainDir$angle
seg1.sel.sp@data$MnDirInv <- seg1.sel.MainDir$angle.inv

seg1.sel.sp@data$MnFl_D <- abs(seg1.sel.sp@data$MnDir - seg1.sel.sp@data$Flow)
seg1.sel.sp@data$MnInfFl_D <- abs(seg1.sel.sp@data$MnDirInv - seg1.sel.sp@data$Flow)
percentage.FlowMinPer <- 20
seg1.sel.sp@data$FlMnPer<- ifelse(((seg1.sel.sp@data$MnFl_D >= (90*((100-percentage.FlowMinPer)/100)) & seg1.sel.sp@data$MnFl_D <= (90*((100+percentage.FlowMinPer)/100)) ) | (seg1.sel.sp@data$MnInfFl_D >= (90*((100-percentage.FlowMinPer)/100)) & seg1.sel.sp@data$MnInfFl_D <= (90*((100+percentage.FlowMinPer)/100)))), 1, 0)


# calculation of length-width-ratio
seg1.sel.LenWidthRatio <- LengthWidthRatio(seg1.sel.sp)
seg1.sel.sp@data$LeWiRat <- seg1.sel.LenWidthRatio$ratio


# transform NA to -9999
replaceInvalids(seg1.sel.sp@data, replace.value = -9999)


# write data to shapefile
rgdal::writeOGR(seg1.sel.sp, dsn = paste(getwd(), dirname(seg1.sel), sep = "/"), overwrite_layer = TRUE,
                layer = file_path_sans_ext(basename(seg1.sel)), driver = "ESRI Shapefile")

# colnames(suppressMessages(shapefiles::read.dbf(paste0(file_path_sans_ext(seg1.sel), ".dbf")))$dbf)
correctDBF(seg1.sel, adjust.n = 25, new.colnames = c("Count", "Slp", "DTM", "DTM_D_51", "Fl_Cos", "Fl_Sin", "SfRghn", "t_Entr", "t_SA", "SVF",
                                                     "Area", "P","P_A", "P_sqrt_A", "Mx_Dist", "D_A", "D_sqrt_A", "Sh_Ind",
                                                     "Flow", "MnDir", "MnDirInv", "MnFl_D", "MnInfFl_D", "FlMnPer", "LeWiRat"))




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ... L1: k-mean thresholding ----------------------------------------
print("... k-mean Thresholding")

# re-read polygons
# rm(seg1.sel.sf)
seg1.sel.sf <- sf::st_read(dsn = paste(getwd(), dirname(seg1.sel), sep = "/"), 
                           layer = file_path_sans_ext(basename(seg1.sel)), quiet = TRUE)


# calculate thresholds
# entropy
l1.thrs.entr <- kMeanThresholds(seg1.sel.sf$t_Entr, size.data = 10)
# l1.thrs.entr.sel <- l1.thrs.entr[2] # empirical chosen
# l1.thrs.entr[order(l1.thrs.entr),]


# dtm to dtm51 dif
# l1.thrs.dtm51dif <- kMeanThresholds(seg1.sel.sf[seg1.sel.sf$t_Entr >= l1.thrs.entr.sel,]$DTM51_Dif)
l1.thrs.dtm51dif <- kMeanThresholds(seg1.sel.sf$DTM_D_51, size.data = 10)
# l1.thrs.dtm51dif.sel <- l1.thrs.dtm51dif[2] # empirical chosen
# l1.thrs.dtm51dif[order(l1.thrs.dtm51dif),]


# sum average
# l1.thrs.sa <- kMeanThresholds(seg1.sel.sf$t_SA)
# l1.thrs.sa <-  kMeanThresholds(seg1.sel.sf[(seg1.sel.sf$t_Entr >= l1.thrs.entr.sel) & (seg1.sel.sf$DTM51_Dif >= l1.thrs.dtm51dif.sel),]$t_SA)
l1.thrs.sa <-  kMeanThresholds(seg1.sel.sf$t_SA, size.data = 10)
# l1.thrs.sa.sel <- l1.thrs.sa[3] # empirical chosen
# l1.thrs.sa[order(l1.thrs.sa),]



# sky view factor
l1.thrs.svf <- kMeanThresholds(seg1.sel.sf$SVF, size.data = 10)
# l1.thrs.svf.sel <- l1.thrs.svf[2] # empirical chosen
# l1.thrs.svf[order(l1.thrs.svf),]


# slope
l1.thrs.slp <- kMeanThresholds(seg1.sel.sf$Slp, size.data = 10)
# l1.thrs.slp.sel <- l1.thrs.slp[3] # empirical chosen
# l1.thrs.slp[order(l1.thrs.slp),]


# CurvMax
l1.thrs.curv15 <- kMeanThresholds(seg1.sel.sf$CurvM15, size.data = 10)


# surface roughness
l1.thrs.srough <- kMeanThresholds(seg1.sel.sf$SfRghn, size.data = 10)
# l1.thrs.srough.sel <- l1.thrs.srough[3] # empirical chosen
# l1.thrs.srough[order(l1.thrs.srough),]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# BREAK ---------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Save Data ---------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
print("... Save Data")
save.image("L1_part1.RData")


