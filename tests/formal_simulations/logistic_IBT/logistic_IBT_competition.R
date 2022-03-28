source("./R/IBTpermutation.R")

betas_size <- c(0, 1, 3)
betas_dist <- c(0, -1, -3)
betas_env <- c(0, .5, 1, 3)

all_betas <- expand.grid(betas_size, betas_dist, betas_env)
all_betas <- as.matrix( cbind(-1, all_betas))

n_rep <- 500
n_res <- length(betas_size) * length(betas_dist) * length(betas_env) * n_rep

res <- data.frame(matrix(NA, n_res, 5))
colnames(res) <- c("p", "beta_size", "beta_dist","beta_env", "i_rep")

set.seed(42)
n_spp <- 20
spp_list <- 1:n_spp

n_island <- 30
island_feature <- data.frame(name = 1:n_island,
                             size = runif(n_island, -1,1),
                             dist = runif(n_island, -1,1))

island_feature2 <- island_feature
island_feature2[,4] <- rnorm(n_island, 0,.5)

spp_trait <- matrix(runif(5*n_spp, -1,1), nrow = n_spp)
trait_dist <- as.matrix(dist(spp_trait))
spp_trait <- data.frame(spp_list, spp_trait)

i_save <- 1
for(i_param in 1:nrow(all_betas)){
  this_beta <- all_betas[i_param,]
  for(i in 1:n_rep){
    cat("setting number: ",i_param, " repeat: ", i, "\n" )
    spp_mat <- simple_competingIBT(spp_list, as.matrix(dist(spp_trait[,2])), island_feature2, this_beta)
    perm_sample <- rho_permutation(spp_mat, trait_dist)
    res[i_save,1] <- 1-perm_sample$p$lower
    res[i_save,2:4] <- this_beta[-1]
    res[i_save,5] <- i
    i_save <- i_save + 1
  }
  write.csv(res, "./tests/formal_simulations/logistic_IBT/p_vals2_competing.csv",row.names = F)
}
