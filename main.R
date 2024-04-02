library(data.table)
library(factoextra)
library(mltools)
library(spatstat)
library(stringr)
library(maptools)
library(virtualspecies)
library(parallel)
library(sf)

source("src/processData.R")
source("src/regression.R")
source("src/maxent.R")
source("src/grid.R")

# projeçoes
UTM22 <- CRS("+proj=utm +zone=22 +south +datum=WGS84 +units=m +no_defs +type=crs")
WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# para reproducibilidade
set.seed(123)

# limites RS e SC
# para uso na geraçao de pontos aleatórios
borders <- shapefile("data/shp/borders/boundaries.shp")

# butiazais
butiazais <- read.csv("data/butia/butia_ponto.csv")
coordinates(butiazais) <- ~x+y
proj4string(butiazais) <- proj4string(borders)

# rasters bioclimáticos e de terreno
rasters <- stack(c(list.files("data/rasters/bioclim", pattern=".tif$", full.names=TRUE),
                   list.files("data/rasters/terrain", pattern=".tif$", full.names=TRUE),
                   list.files("data/rasters/site_dist", full.names=T, pattern=".tif$")))

# remoçao de colinearidade usando um limite de 0.7
removeCollinearity(rasters, multicollinearity.cutoff=0.7, sample.points=TRUE,
                   nb.points=100000, plot=TRUE)
dev.print(pdf, "figs/cluster.pdf")
dev.off()

SELECTED_VARS <- c("bio1", "bio10", "bio12", "bio13", "bio15", "bio2",
                   "bio3", "bio4", "bio8", "bio9", "dcoast", "drivs", "slo",
                   "aceramic_m_dist", "all_sites_m_dist", "cerritos_m_dist",
                   "guarani_m_dist")

SELECTED_VARS_ENV <- c("bio1", "bio10", "bio12", "bio13", "bio15", "bio2",
                       "bio3", "bio4", "bio8", "bio9", "dcoast", "drivs", "slo")

# retemos apenas as camadas selecionadas
all_layers <- rasters[[which(names(rasters) %in% SELECTED_VARS)]]
plot(all_layers)
dev.print(pdf, "figs/all_layers.pdf")
dev.off()

layers_env <- rasters[[which(names(rasters) %in% SELECTED_VARS_ENV)]]

# conversao das variáveis ambientais a UTM22
covars <- projectRaster(layers_env, res=10000, crs=UTM22)

# convertemos de raster a image pixels e criamos uma lista
gvars <- list()
for (i in 1:nlayers(covars)) {
    gvars[[i]] <- as.im(as(covars[[i]], "SpatialGridDataFrame"))
    names(gvars)[[i]] <- names(covars)[i]
}

# fórmula com as variáveis independentes: ~wc2.1_30s_bio_12 + wc2.1_30s_bio_15 + wc2.1_30s_bio_18 + ...
gtrend <- reformulate(names(layers_env))


# fill.na <- function(x, i=13) {
#     if (is.na(x)[i])
#         return (mean(x, na.rm=TRUE))
#     else
#         return (x[i])
# }


# BLR
points <- make_rnd_data(butiazais, borders, num_rand=500)

vals <- extract_env_data(points, all_layers)
k.nn <- spdep::knn2nb(spdep::knearneigh(vals@coords, longlat=T))
bw <- max(unlist(spdep::nbdists(k.nn, vals@coords)))
vals$AutoCov <- spdep::autocov_dist(vals$presence, xy = vals@coords, nbs = bw, style = "W", type = "inverse")

blr_model <- binary_train(vals@data)
sink("tables/blr.txt")
print(summary(blr_model$model))
print(blr_model$cv_score)
sink()

# MaxEnt

dir.create("tables/maxent")

max_model <- maxent_train(points, all_layers)   # AUC de ~93%
plot(max_model$model)
dev.print(pdf, "figs/maxent.pdf")
dev.off()

#plot maxent prediction
pred <- predict(max_model$model, all_layers)
plot(pred)
plot(butiazais, add=T)
dev.print(pdf, "figs/maxent_pred.pdf")
dev.off()

run_ripley_tests <- function(site_name="all_sites") {
    # convertemos os pontos dos butiazais para UTM e definimos uma coluna
    # com o tipo (marca) dos pontos
    butias_m <- st_transform(st_as_sf(butiazais), UTM22)
    butias_m$type <- "butia"
    butias_m$type <- as.factor(butias_m$type)

    # mesmo procedimento para os sítios arqueológicos
    # o arquivo contém uma filtragem dos sítios com no máximo 1 sítio por km^2
    # (aproximadamente a resoluçao dos rasters)
    sites <- read.csv(paste("data/sites/filtrado/", site_name, ".csv", sep=""), header = TRUE, sep=",")
    coordinates(sites) <- ~decimallon+decimallat
    proj4string(sites) <- WGS84
    sites_m <- st_transform(st_as_sf(sites), UTM22)
    sites_m$type <- "site"
    sites_m$type <- as.factor(sites_m$type)
   
    # convertemos para um objeto point pattern
    df <- data.frame(x=c(st_coordinates(butias_m)[,1], st_coordinates(sites_m)[,1]), y=c(st_coordinates(butias_m)[,2], st_coordinates(sites_m)[,2]),
                    type=c(butias_m$type, sites_m$type))
    marked <- ppp(df$x, df$y, marks=df$type, window=ripras(df$x,df$y))

    #plot(marked)

    # modelo básico - sem covariáveis ambientais
    basic <- envelope(marked, fun=pcfcross, nsim=999, correction="iso", divisor="d")
    plot(basic)
    dev.print(pdf, paste("figs/ripley_basic_", site_name, ".pdf", sep=""))
    dev.off()

    fit <- step(ppm(marked, trend=gtrend, covariates=gvars, correction="iso"),
                direction="both", k=2)
    sink(paste("tables/ppm_", site_name, ".txt", sep=""))
    print(summary(fit))
    sink()

    # modelo de "primeira ordem" - fatóres "exógenos" - ao invés de uma distribuiçao
    # aleatório, comparamos com a distribuiçao prevista a partir das covariáveis
    # ambientais
    simfit <- envelope(fit, fun=pcfcross, nsim=999, correction="iso", divisor="d")
    plot(simfit)
    dev.print(pdf, paste("figs/ripley_env_", site_name, ".pdf", sep=""))
    dev.off()
}

site_names <- c("all_sites", "guarani", "proto_je", "aceramic", "cerritos", "sambaquis")
for (site_name in site_names) {
  run_ripley_tests(site_name)
}
