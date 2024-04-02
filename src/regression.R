library(MASS)
library(glmnet)
library(spatialEco)
library(caret)
library(pROC)

# ---------------------------------------------------------
# Funçoes para executar a regressao binária logística (BLR)
# ---------------------------------------------------------

# Executa regressa binária logística
#
# Argumentos
# ----------
#   x : pontos de presença e ausência; data frame com, minimamente, uma coluna
#       "presence" com valores 1 (presença) e 0 ([pseudo-]ausência) e colunas
#       para variáveis independentes.
#   test_size : percentual da amostra a ser separada para teste
#
binary_train <- function(x) {
    folds <- createFolds(x$presence, 5)
    scores <- c()
    for (i in 1:5)
    {
        X_train <- x[-folds[[i]],]
        X_test <- x[folds[[i]],]

        # modelo com todas as variáveis
        full_model <- glm(presence ~ ., family=binomial(link="logit"), data=X_train)

        # variáveis com pouca significância ou colinearidade sao removidas até
        # se chegar ao modelo mais parsimonioso possível e com mais alto AIC
        step_model <- stepAIC(full_model, direction="both", trace=FALSE)

        y_hat <- predict(step_model, X_test, type="response")
        score <- roc(response=X_test$presence, predictor=y_hat)
        scores[i] <- score$auc
    }

    final_model_full <- glm(presence ~ ., family=binomial(link="logit"), data=x)
    final_model_step <- stepAIC(final_model_full, direction="both", trace=FALSE)

    print(scores)
    return(list(cv_score=mean(scores), model=final_model_step))
}
