
qdh_traits <- apply(qdh_traits, 2, function(w){(w-mean(w))/sd(w)})
trait_dist1 <- as.matrix(dist(qdh_traits[,1])) # body size
trait_dist2 <- as.matrix(dist(qdh_traits[,-1])) # other than body size

par(mfrow = c(1,2))
perm_sample1 <- rho_permutation(qdh_spp_mat, trait_dist1, n = 3000)
plot(density(perm_sample1$permutations))
abline(v = perm_sample1$observed, col = "red")
perm_sample1 <- rho_permutation(qdh_spp_mat, trait_dist2, n = 3000)
plot(density(perm_sample2$permutations))
abline(v = perm_sample2$observed, col = "red")

perm_sample1$p
perm_sample2$p
