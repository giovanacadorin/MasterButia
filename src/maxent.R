library(dismo)

# ---------------------------------------------------------
# Funçoes para executar o modelo de máxima entropia
# ---------------------------------------------------------

# Executa MaxEnt
#
# Argumentos
# ----------
#   x : um objeto em formato SpatialPoints* com uma coluna binária presence
#       onde 1 sao as observaçoes
#   raster_stack : um stack de rasters com as covariáveis ambientais
#   factors : nome das variáveis categóricas no stack de rasters
#   test_size : percentual de pontos a serem separados para teste
#
maxent_train <- function(x, raster_stack, factors=NULL, test_size=0.2) {
    x_pres <- x[x$presence==1,]
    x_abs <- x[x$presence==0,]
    maxent_model <- maxent(raster_stack, x_pres, args=c("responsecurves=true"),
                           factors=factors, path="tables/maxent")
    eval <- evaluate(x_pres, x_abs, maxent_model, raster_stack)
    print(paste("AUC:", eval@auc))
    return(list(model=maxent_model, eval=eval))
}
