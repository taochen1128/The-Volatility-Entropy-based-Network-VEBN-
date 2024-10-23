returns_data<-read.csv("return.csv", header = TRUE, sep = ",")

##garch(1.1)--------------------------------
#Define the GARCH(1,1) model specification
spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                   mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                   distribution.model = "std")

#Fit the GARCH model to each time series
garch_fit_1 <- ugarchfit(spec = spec, data = returns_data[,1])
garch_fit_2 <- ugarchfit(spec = spec, data = returns_data[,2])
garch_fit_3 <- ugarchfit(spec = spec, data = returns_data[,3])

# Extracted conditional volatility
volatility_1 <- sigma(garch_fit_1)
volatility_2 <- sigma(garch_fit_2)
volatility_3 <- sigma(garch_fit_3)
volatility_df <- data.frame(
  Volatility_1 = volatility_1,
  Volatility_2 = volatility_2,
  Volatility_3 = volatility_3,)

##bayesian Stochastic Volatility Estimates------------------------------------
volatility_list <- list()

for (i in 1:ncol(returns_data)) {
  selected_returns <- returns_data[, i]
  selected_returns <- as.numeric(selected_returns)
  selected_returns <- selected_returns[!is.na(selected_returns)]
  res <- svsample(selected_returns, priormu = c(-10, 1), priorphi = c(20, 1.5), priorsigma = 0.1)
  sd_exp_ht <- as.numeric(summary(res, showlatent = TRUE)$latent[, 7])
  volatility_list[[colnames(returns_data)[i]]] <- sd_exp_ht
}

print(volatility_list)

##garch(p,q)------------------------------------
GARCH.likelihood.Forecast <- function(data, max_lag = 20, significance_level = 0.05) {
  auto_select_garch_order <- function(data, max_lag, significance_level) {
    acf_values <- acf(data, plot = FALSE)$acf[-1]
    pacf_values <- pacf(data, plot = FALSE)$acf
    n <- length(data)
    threshold <- qnorm(1 - significance_level / 2) / sqrt(n)
    significant_p <- which(abs(pacf_values) > threshold)
    significant_q <- which(abs(acf_values) > threshold)
    max_p <- ifelse(length(significant_p) > 0, min(max(significant_p), max_lag), 1)
    max_q <- ifelse(length(significant_q) > 0, min(max(significant_q), max_lag), 1)
    return(list(p = max_p, q = max_q))
  }
  
  garch_volatility_list <- list()
  garch_order_list <- list()
  
  for (i in 1:ncol(data)) {
    datasetv <- as.ts(data[, i])
    pqvalue <- auto_select_garch_order(datasetv, max_lag, significance_level)
    p <- as.numeric(pqvalue$p)
    q <- as.numeric(pqvalue$q)
    
    spec <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(p, q)),
      mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
      distribution.model = "norm"
    )
    
    garch_fit <- ugarchfit(spec = spec, data = datasetv)
    
    
    fitted_volatility <- sigma(garch_fit)
    
    garch_volatility_list[[colnames(data)[i]]] <- fitted_volatility
    garch_order_list[[colnames(data)[i]]] <- c(p, q)
  }
  
  garch_volatility_df <- do.call(cbind, garch_volatility_list)
  garch_order_df <- do.call(rbind, garch_order_list)
  colnames(garch_order_df) <- c("p", "q")
  
  return(list(volatility = garch_volatility_df, order = garch_order_df))
}

result <- GARCH.likelihood.Forecast(returns_data)


#Calculate the effective transfer entropy-----------------------------
tetransferen <- function(d) {
  
  te <- transfer_entropy(d[,1], d[,2], shuffles = 50, nboot = 100, quiet = T)
  coef(te)[1:2, 2:3]
  
}

ete.result <- data.frame(para = matrix(0,nrow = length(volatility_df),ncol = length(volatility_df)))  
names(ete.result) <- names(volatility_df)
row.names(ete.result) <- names(volatility_df)

for (i in 1:(length(volatility_df)-1)) {
  for (j in (i+1):length(volatility_df)) {
    outv <- tetransferen(volatility_df[,c(i,j)])
    ete.result[i,j] <- outv[1,1]
    ete.result[j,i] <- outv[2,1]  
    #se.tetransferen[i,j] <- outv[1,2]
    #se.tetransferen[j,i] <- outv[2,2]    
  }
}

#paiewise ETE-----
te_vo<-read.csv("ete.result.csv",header=TRUE,sep=",")
te_vo=te_vo[,-1]

pairwise.Ete_vo <- data.frame(para = matrix(0,nrow = length(te_vo),
                                            ncol = length(te_vo)))  
names(pairwise.Ete_vo) <- names(te_vo)
row.names(pairwise.Ete_vo) <- names(te_vo)
###
for (i in 1:(length(te_vo)-1)) {
  for (j in (i+1):length(te_vo)) {
    pairwise.Ete_vo[i,i] <-0
    pairwise.Ete_vo[i,j] <- te_vo[i,j]-te_vo[j,i] 
    pairwise.Ete_vo[j,i] <- te_vo[j,i]-te_vo[i,j]  
    if (pairwise.Ete_vo[i,j] < 0) {
      pairwise.Ete_vo[i,j] <- 0
    } else {
      pairwise.Ete_vo[i,j] <- pairwise.Ete_vo[i,j]
      # Code to record the value in your table goes here
    }
    if (pairwise.Ete_vo[j,i] < 0) {
      pairwise.Ete_vo[j,i] <- 0
    } else {
      pairwise.Ete_vo[j,i] <- pairwise.Ete_vo[j,i]
      # Code to record the value in your table goes here
    }
    
  }
}

# Make sure the data is read correctly as a matrix
matrix_result<- as.matrix(pairwise.Ete_vo)

#Creating the igraph network object
g <- graph_from_adjacency_matrix(matrix_result, mode = "directed", weighted = TRUE)

# Calculating thresholds--------------------------------------------
threshold <- quantile(E(g)$weight, 0.70)#80%-90%-mean

# Creating a threshold network
g_thresholded <- delete_edges(g, E(g)[E(g)$weight < threshold])

# Calculating degree
node_degree <- degree(g_thresholded)

x=node_degree
x <- x[x > 0]

if(length(x) > 0){
  m_pl <- displ$new(x)
} else {
  stop("No data left after removing isolated nodes.")
}

ticks_x = 5
ticks_y = 5 

fit <- fit_power_law(x,force.continuous = FALSE)

#fit <- fit_power_law(x,xmin = 22,force.continuous = FALSE)
#m_pl$setXmin(min(x))     
m_pl$setXmin(fit$xmin)                      
m_pl$setPars(fit$alpha)
xmin <- m_pl$getXmin()
gma <- m_pl$getPars()

cdf <- plot(m_pl, draw = FALSE)

fit[["alpha"]]
fit[["KS.stat"]]
fit[["KS.p"]]

####plot----------------------------------------
lfit <- tibble(k = 1:max(x), p = (1:max(x) / xmin)^(1 - gma))
lfit$p <- lfit$p * cdf$y[which(cdf$x >= xmin)[1]]

mean_pl <- ggplot(cdf, aes_(~x, ~y)) + geom_point() +
  geom_line(data = lfit, aes_(~k, ~p), colour = "red") +
  coord_cartesian(xlim = c(min(cdf$x), max(cdf$x)),
                  ylim = c(min(cdf$y), 1)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = ticks_x), 
                labels = trans_format("log10", math_format())) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = ticks_y),
                labels = trans_format("log10", math_format())) + 
  annotate(
    label = "p-value ==0.8634 ",parse = TRUE,size=4) + labs(x = "IN(k)(degree_70%)", y = "IN(P(X))") +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(family = "Times New Roman", size = 25),
                     axis.title = element_text(family = "Times New Roman", size = 25),
                     plot.title = element_text(family = "Times New Roman", size = 20))
plot(mean_pl)

### Calculating network features------------------------------
# Average path length
average_path_length <- mean_distance(g_thresholded, directed = TRUE)

# clu factor
clustering_coefficient <- transitivity(g_thresholded, type = "average")

#Average degree
average_degree <- mean(degree(g_thresholded))

# Network density
network_density <- graph.density(g_thresholded)

#DC
edge_weights <- E(g_thresholded)$weight
mean_edge_weight <- mean(edge_weights, na.rm = TRUE)
var_edge_weight <- var(edge_weights, na.rm = TRUE)
sd_edge_weight <- sd(edge_weights, na.rm = TRUE)

#RC
node_degrees <- degree(g_thresholded)
mean_node_degree <- mean(node_degrees)
var_node_degree <- var(node_degrees)
sd_node_degree <- sd(node_degrees)

###-nodal centrality-------------------------------------------------
node_degree <- degree(g_thresholded)
node_closeness <- closeness(g_thresholded, mode = "all")
node_betweenness <- betweenness(g_thresholded, directed = TRUE)
node_eigen <- eigen_centrality(g_thresholded, directed = TRUE)$vector
node_pagerank <- page_rank(g_thresholded)$vector

