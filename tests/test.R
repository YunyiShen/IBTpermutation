source("./R/IBTpermutation.R")
n_spp <- 15
spp_list <- 1:n_spp

n_island <- 25
island_feature <- data.frame(name = 1:n_island,
                             size = runif(n_island, -1,1),
                             dist = runif(n_island, -1,1))

spp_trait <- matrix(runif(5*n_spp, -1,1), nrow = n_spp)
trait_dist <- as.matrix(dist(spp_trait))
spp_trait <- data.frame(spp_list, spp_trait)

# IBT
ps_null <- rep(NA, 100)
for(i in 1:100){
  betas <- c(-1,.5,-.5)
  spp_mat <- simple_IBT(spp_list, island_feature, betas)
  perm_sample <- rho_permutation(spp_mat, trait_dist)
  #hist(perm_sample$permutations)
  #abline(v = perm_sample$observed, col = "red")
  ps_null[i] <- perm_sample$p
  #quantile(perm_sample$permutations, c(.025,.975))
}
hist(ps_null)


# IBT with sorting
ps_h1 <- rep(NA, 100)
island_feature2 <- island_feature
island_feature2[,4] <- rnorm(n_island, 0,.5)
for(i in 1:100){
  betas_sorting <- c(-1,0,-0, -1)
  spp_mat_sorting <- simple_sortingIBT(spp_trait, island_feature2, betas_sorting)
  perm_sample <- rho_permutation(spp_mat_sorting, trait_dist, n = 3000)
  hist(perm_sample$permutations)
  abline(v = perm_sample$observed, col = "red")
  ps_h1[i] <- perm_sample$p
  quantile(perm_sample$permutations, c(.025,.975))
}

plot(ecdf(ps_null))
plot(ecdf(ps_h1))

