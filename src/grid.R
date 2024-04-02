library(exactextractr)
library(nlme)

# ---------------------------------------------------------
# Funçoes para executar a regressao GLS a partir dos valores
# médios em células de uma grade n x n km
# ---------------------------------------------------------

# Extrai os valores da grade para um data frame.
#
# Argumentos
# ----------
#   grid : uma grade em formato shapefile
#   layers : um stack de rasters ambientais com os valores a serem extraídos
#
get_grid_data <- function(grid, layers) {
    values <- exact_extract(layers, grid, fun="mean", progress=FALSE)
    grid@data <- cbind(grid@data, values)
    grid$x <- grid@data$left
    grid$y <- grid@data$top
    gc()
    # as primeiras 5 colunas sao eliminadas (id, left, top, etc.)
    grid@data <- grid@data[, -(1:5)]
    return (grid@data)
}

# Executa regressao GLS
#
# Argumentos
# ----------
#   grid : uma grade em formato shapefile com as colunas num_butia para número
#          de butiazais (variável dependente) e num_sites para número de sítios
#          arqueológicos por célula
#   layers : um stack de rasters ambientais com as variáveis independentes
#
run_gls_model <- function(grid, layers) {
    data <- get_grid_data(grid, layers)
    # um modelo espacial é incluído através do parâmetro correlation para
    # dar conta da autocorrelaçao espacial
    model <- gls(log(num_butia+1)~.-x-y, data=data, na.action=na.omit,
                 correlation=corSpher(form=~x+y, nugget=TRUE))
    print(summary(model)$tTable)
    return(list(data=data, model=model))
}
