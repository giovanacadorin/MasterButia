library(raster)

# ---------------------------------------------------------
# Funçoes gerais de preparaçao dos dados
# ---------------------------------------------------------

# Cria pontos aleatórios de pseudo-ausência
#
# Argumentos
# ----------
#   points : um objeto de tipo SpatialPoints* com os pontos de presença
#   borders : shapefile com os limites da área de estudo
#   num_rand : número de pontos aleatórios a serem gerados (opcional)
#
make_rnd_data <- function(points, borders, num_rand=NULL) {
    # se nao é especificado o número de pontos, geram-se tantos quantos forem
    # os pontos de presença
    if (is.null(num_rand)) size <- length(points)
    else size <- num_rand

    rnd <- as(spsample(borders, size, "random", seed=123), "SpatialPointsDataFrame")

    # uma coluna presence é criada com valores 0 para os pontos gerados de pseudo-
    # ausência e 1 para os pontos de presença
    rnd@data <- data.frame(presence=rep(0, length(rnd)))
    points$presence <- 1

    return(bind(points, rnd))
}

# Extrai valores de um stack de rasters para os pontos gerados na funçao anterior
#
# Argumentos
# ----------
#   points : pontos de presença e ausência (SpatialPointsDataFrame)
#   rasters : um stack de rasters com os valores a serem extraídos
#
extract_env_data <- function(points, rasters) {
    values <- extract(rasters, points)
    df <- as.data.frame(cbind(values, points$presence))
    colnames(df)[ncol(df)] <- "presence"
    coordinates(df) <- coordinates(points)
    df <- df[!is.na(df@data[,1]),]
    return(df)
}