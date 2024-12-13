# Libraries ----
library(dplyr)
library(reshape2)
library(Matrix)
library(glmnet)
library(stringr)
library(igraph)

source("./Code/functions.R")
devtools::load_all("./Code/bsglm/")

# Build (league) Design Matrix ----
file <- "./Data/NBA-PbP-Sample-Dataset.csv"

df <- read.csv(file)

# only take columns we need and rename for ease
df <- df[2:nrow(df)
         , c(1       # game ID
             , 4:13  # team lineup
             , 15:16 # home and away aggregate scores
             , 19    # duration
             , 33)]  # points scored on play
names(df) <- c("game", paste0("a", 1:5), paste0("h", 1:5), "away", "home", "time", "change")
df$change <- as.numeric(df$change)
df$change[is.na(df$change)] <- 0
df$home <- as.numeric(df$home)
df$away <- as.numeric(df$away)
df$lead <- df$home - df$away
df$time <- as.numeric(substr(df$time, 6, 7))

# track lead change
df$previous_lead <- c(0, df$lead[-nrow(df)])
df$pm <- df$change * (2*(df$lead > df$previous_lead) - 1)

# fix lex order
df[, 2:6] <- t(apply(df[, 2:6], 1, sort))
df[, 7:11] <- t(apply(df[, 7:11], 1, sort))

# group by lineup
df <- df %>%
  group_by(game, a1, a2, a3, a4, a5, h1, h2, h3, h4, h5) %>%
  dplyr::summarise(pm = sum(pm)
                   , time = sum(time))

df$lineup <- rownames(df)

players_dal <- unique(unlist(df[, 2:6]))
players_bos <- unique(unlist(df[, 7:11]))

# make design matrix
X <- melt(df, id.vars = c("lineup", "pm", "time", "game"))
X <- X[, c(1, 2, 3, 6)]

Y <- unique(X[, c(1, 2)])
W <- unique(X[, c(1, 3)])

X <- table(X[, c(1, 4)]) > 0

# match X, Y, and W
Y <- Y$pm[match(as.numeric(rownames(X)), Y$lineup)]
W <- W$time[match(as.numeric(rownames(X)), W$lineup)]

# remove 0 second lineups
X <- X[W > 0, ]
Y <- Y[W > 0]
W <- W[W > 0]

X <- X[, colSums(abs(X)) > 0] # substitutions can have blank players

players <- colnames(X)
away_players <- unique(unlist(df[, paste0("a", 1:5)]))
away_ind <- players %in% away_players
X[, players %in% away_players] <- -X[, away_ind]

X1 <- model.matrix(~.^2, data = as.data.frame(X))[,-1]
all_pairs <- apply(combn(players, 2), 2, paste, collapse = ":")
colnames(X1) <- c(players, all_pairs)

# remove pairs between teams
team_pairs <- c(apply(combn(players[away_ind], 2), 2, paste, collapse = ":")
                , apply(combn(players[!away_ind], 2), 2, paste, collapse = ":"))
X1 <- X1[, which(colnames(X1) %in% c(players, team_pairs))]

# remove pairs that never played together
X1 <- X1[, colSums(abs(X1)) > 0]

# away pairs need -1
for (j in (length(players)+1):ncol(X1)) {
  p1 <- gsub(":.*", "", colnames(X1)[j])
  p2 <- gsub(".*:", "", colnames(X1)[j])
  X1[, j] <- X1[, j] * pmin(X1[, p1], X1[, p2])
}

X <- Matrix(as.matrix(X), sparse = TRUE)
X1 <- Matrix(as.matrix(X1), sparse = TRUE)

# PM ----
PM <- Y %*% X

# APM_league ----
APM_RR <- cv.glmnet(X, Y, weights = W, alpha = 0, intercept = FALSE)
APM_RR <- APM_RR$glmnet.fit$beta[, which(APM_RR$lambda == APM_RR$lambda.min)]

# PAPM_league ----
PAPM_RR <- cv.glmnet(X1, Y, weights = W, alpha = 0, intercept = FALSE)

df_dal <- cbind(PM[match(players_dal, colnames(PM))]
                , APM_RR[players_dal])

df_bos <- cbind(PM[match(players_bos, colnames(PM))]
                , APM_RR[players_bos])

# By Team ----
for (team in 0:1) {
  df <- read.csv(file)
  
  # only take columns we need and rename for ease
  df <- df[2:nrow(df)
           , c(1                      # game ID
               , 4:8 + 5*team         # team lineup
               , 15:16                # home and away aggregate scores
               , 19                   # duration
               , 33)]                 # points scored on play
  names(df) <- c("game", paste0("p", 1:5), "away", "home", "time", "change")
  df$change <- as.numeric(df$change)
  df$change[is.na(df$change)] <- 0
  df$home <- as.numeric(df$home)
  df$away <- as.numeric(df$away)
  df$lead <- (df$home - df$away) * (2*team - 1)
  df$time <- as.numeric(substr(df$time, 6, 7))
  
  # track lead change
  df$previous_lead <- c(0, df$lead[-nrow(df)])
  df$pm <- df$change * (2*(df$lead > df$previous_lead) - 1)
  
  # fix lex order
  df[, 2:6] <- t(apply(df[, 2:6], 1, sort))
  
  # group by lineup
  df <- df %>%
    group_by(game, p1, p2, p3, p4, p5) %>%
    dplyr::summarise(pm = sum(pm)
                     , time = sum(time))
  
  df$lineup <- rownames(df)
  
  # Build (team) Design Matrix ----
  # make design matrix
  X <- melt(df, id.vars = c("lineup", "pm", "time", "game"))
  X <- X[, c(1, 2, 3, 6)]
  
  Y <- unique(X[, c(1, 2)])
  W <- unique(X[, c(1, 3)])
  
  X <- table(X[, c(1, 4)]) > 0
  
  # match X, Y, and W
  Y <- Y$pm[match(as.numeric(rownames(X)), Y$lineup)]
  W <- W$time[match(as.numeric(rownames(X)), W$lineup)]
  
  # remove 0 second lineups
  X <- X[W > 0, ]
  Y <- Y[W > 0]
  W <- W[W > 0]
  
  # PM ----
  PM <- Y %*% X
  
  # APM ----
  APM_RR <- cv.glmnet(X, Y, weights = W, alpha = 0, intercept = FALSE)
  APM_RR <- APM_RR$glmnet.fit$beta[, which(APM_RR$lambda == APM_RR$lambda.min)]
  
  # Build (pairs) Design Matrix ----
  X1 <- model.matrix(~.^2, data= as.data.frame(X))[,-1]
  pairs <- apply(combn(colnames(X), 2), 2, paste, collapse = ":")
  colnames(X1) <- c(colnames(X), pairs)
  
  # remove pairs that never played together
  X1 <- X1[, colSums(X1) > 0]
  
  # pairs APM
  PAPM_RR <- cv.glmnet(X1, Y, weights = W, alpha = 0, intercept = FALSE)
  PAPM_RR <- PAPM_RR$glmnet.fit$beta[, which(PAPM_RR$lambda == PAPM_RR$lambda.min)]
  
  # HAPM ----
  players <- colnames(X)
  combns <- c(players
              , apply(combn(players, 2), 2, paste, collapse = ":") #could just  use pairs, but lets assume I dont run that
              , apply(combn(players, 3), 2, paste, collapse = ":")
              , apply(combn(players, 4), 2, paste, collapse = ":")
              , apply(combn(players, 5), 2, paste, collapse = ":"))
  X2 <- matrix(0, nrow = length(combns), ncol = length(players))
  Y2 <- W2 <- vector(length = length(combns))
  colnames(X2) <- players
  rownames(X2) <- names(Y2) <- names(W2) <- combns
  for (i in seq_len(length(combns))) {
    combo <- unlist(strsplit(combns[i], ":"))
    player.ind <- which(players %in% combo)
    X2[i, player.ind] <- 1
    Y2[i] <- as.numeric(apply(X[, player.ind, drop = FALSE], 1, prod) %*% Y)
    W2[i] <- as.numeric(apply(X[, player.ind, drop = FALSE], 1, prod) %*% W)
  }
  
  # remove combinations that never played together
  combns <- combns[W2 > 0]
  X2 <- X2[W2 > 0, ]
  Y2 <- Y2[W2 > 0]
  W2 <- W2[W2 > 0]
  
  # hypergraph APM
  HAPM_RR <- cv.glmnet(X2, Y2, weights = W2, alpha = 0, intercept = FALSE)
  
  # LAPM ----
  # lets make the adjacency matrix
  A <- as.matrix(proxy::dist(str_split(rownames(X2), ":"), method = jaccard_dist))
  
  # Now that you have the adjacency matrix - create the graph/laplacian/eigenvectors
  g <- graph_from_adjacency_matrix(A, weighted = TRUE, mode = "undirected")
  
  # Laplacian of l
  lap_m <- laplacian_matrix(g)
  eig <- eigen(lap_m)
  p <- dim(X2)[1] #number of nodes
  design <- eig$vectors[,p:1] #changing order of eigenvectors
  e.values <- rev(eig$values) #changing order of eigenvalues
  
  dim <- elbow_finder(seq(1:length(e.values[-1])), 1/e.values[-1]) 
  expand <- dim[1] + 1
  design <- design[,2:expand] # remove intercept, because bglm.fit adds it back in. 
  
  lam <- 0.0001 # given an expansion, this doesn't matter much (because we care about ranks)
  
  ##update lines below to match lines above
  prior.coef <- list(mean=rep(0, length = expand), precision=lam*e.values[1:expand])
  LAPM_all <- glm(Y2 ~ design, weights = W2, family = "gaussian", method = bsglm::fitter(prior_coef = prior.coef, adjust_dispersion = TRUE))
  pred <- LAPM_all$fitted.values
  LAPM_pred <- data.frame(Player = rownames(X2), LAPM_team = pred)
  
  if (team) {
    df_bos <- cbind(df_bos
                    # , PM[match(players_bos, colnames(PM))]
                    , APM_RR[players_bos]
                    , HAPM_RR$glmnet.fit$beta[!grepl(":", colnames(X2)), which(HAPM_RR$lambda == HAPM_RR$lambda.min)][players_bos]
                    , LAPM_pred[players_bos, 2]
    )
  } else {
    df_dal <- cbind(df_dal
                    # , PM[match(players_dal, colnames(PM))]
                    , APM_RR[players_dal]
                    , HAPM_RR$glmnet.fit$beta[!grepl(":", colnames(X2)), which(HAPM_RR$lambda == HAPM_RR$lambda.min)][players_dal]
                    , LAPM_pred[players_dal, 2]
    )
  }
}

colnames(df_bos) <- colnames(df_dal) <- c("PM"
                                          , "APM_league"
                                          , "APM"
                                          , "HAPM"
                                          , "LAPM")
