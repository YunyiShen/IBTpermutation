source("./R/IBTpermutation.R")
#qdh_traits <- apply(qdh_traits, 2, function(w){(w-mean(w))/sd(w)})
trait_dist1 <- as.matrix(dist(qdh_traits[,1])) # body size
trait_dist2 <-  0 * trait_dist1 # other than body size

# calculate a simple niche distance, i.e. if there are overlap then it is 0 otherwise it is 1.
for(i in 1:nrow(trait_dist2)){
  for(j in 1:ncol(trait_dist2)){
    temp <- qdh_traits[c(i,j),-1]
    trait_dist2[i,j] <- 2-max(colSums(temp))
  }
}


set.seed(12345)
par(mfrow = c(1,2))
perm_sample1 <- rho_permutation(qdh_spp_mat, trait_dist1, n = 3000)
plot(density(perm_sample1$permutations), main = "body size", xlab = "beta")
abline(v = perm_sample1$observed, col = "red")
perm_sample2 <- rho_permutation(qdh_spp_mat, trait_dist2, n = 3000)
plot(density(perm_sample2$permutations), main = "niche", xlab = "beta")
abline(v = perm_sample2$observed, col = "red")

perm_sample1$p
perm_sample2$p
