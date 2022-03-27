sigmoid <- function(x){
  1/(1+exp(-x))
}

#' simulating a simple IBT presence matrix
#' simulating a species presence matrix using logistic presence probability depends on distance and size of island
#' @param spp_list string vector of species
#' @param island_feature a matrix with three columns, name of island, size of island, and distance to mainland
#' @param betas the coefficients for the logistics presence probability, intercept, size and distance
#' @return a species presence matrix simulated using IBT, row as islands
simple_IBT <- function(spp_list, island_feature, betas){
  designmat <- apply(as.matrix(island_feature[,2:3]),2,
                     function(w){(w-mean(w))/sd(w)})
  ps <- cbind(1, designmat) %*% betas |> sigmoid()

  presence_mat <- matrix(runif(length(spp_list)*nrow(island_feature)), ncol = length(spp_list)) |>
    apply(2,function(w,p){w<=p},ps)
  rownames(presence_mat) <- island_feature[,1]
  colnames(presence_mat) <- spp_list
  return(presence_mat)

}



#' simulating a simple IBT presence matrix with simple environment sorting
#' simulating a species presence matrix using logistic presence probability depends on distance and size of island, while a island feature and a species trait, when they are close the species are more likely to be present
#' @param spp_traits trait of species, should be a matrix, first column should be name and second being the trait being selected
#' @param island_feature a matrix with at least four columns, name of island, size of island, and distance to mainland, the last one is the target trait of an island, the distance between the required trait and the trait of the species determined the chance (together with size and distance)
#' @param betas the coefficients for the logistics presence probability, intercept, size, distance and island requirement difference to trait of species
#' @return a species presence matrix simulated using IBT, row as islands
simple_sortingIBT <- function(spp_traits, island_feature, betas){
  designmat <- apply(as.matrix(island_feature[,2:3]),2,
                     function(w){(w-mean(w))/sd(w)})

  sorting_mat <- abs(island_feature[,4] %*% t(rep(1, nrow(spp_traits))) - rep(1, nrow(island_feature)) %*% t(spp_traits[,2]))

  sorting_mat <- (sorting_mat-mean(sorting_mat))/sd(sorting_mat)
  ps <- (cbind(1, designmat) %*% betas[-4] %*% t(rep(1, nrow(spp_traits))) +
    betas[4] *  sorting_mat) |>
    sigmoid()
  #browser()
  presence_mat <- matrix(runif(nrow(spp_traits)*nrow(island_feature)), ncol = nrow(spp_traits)) <= ps

  rownames(presence_mat) <- island_feature[,1]
  colnames(presence_mat) <- spp_traits[,1]
  return(presence_mat)

}

#' test for non-neutral factor using rho statistic and its permutation distribution
#' rho statistics is the least square coefficient between richness and average trait distance change on an island from grand mean distance (when the island has more than one species), this routine also draw a permutation distribution by switching species labels to test for competition (more positive) and environmental sorting (more negative)
#' @param spp_mat presence matrix, row as islands/sites, column as species
#' @param trait_dist trait distance matrix between species, order should be the same as spp_mat's column
#' @param n number of random permutations to draw
#' @return a list with the observed rho, permutation samples and one side/double side p-value
rho_permutation <- function(spp_mat, trait_dist, n = min(2000, gamma(ncol(spp_mat)+1))){
  spp_mat <- spp_mat==1
  permutations <- lapply(1:n, function(i,spp_mat){sample(1:ncol(spp_mat))},spp_mat)
  permutations <- c(list(1:ncol(spp_mat)),permutations)
  temp <- sapply(permutations,function(perm, spp_mat, trait_dist){
    temp <- matrix(0, nrow = nrow(spp_mat), ncol = 2)
    temp[,2] <- rowSums(spp_mat) # richness
    for(i in 1:nrow(spp_mat)){
      if(temp[i,2]<=1) next
      perm_temp <- trait_dist[perm,perm] # get the permuted dist matrix, this is a bad way of implementation as most of the distances are not needed
      perm_temp <- perm_temp[spp_mat[i,],spp_mat[i,]]
      temp[i,1] <- (mean(perm_temp[upper.tri(perm_temp)])) # change in  average distance
    } # loop over each site
    #browser()
    temp[,1] <- (temp[,1] - mean(trait_dist[upper.tri(trait_dist)]))/sd(trait_dist[upper.tri(trait_dist)]) # expected change from grand mean
    #temp[,1] <- (temp[,1] -temp[1,1])
    #return(cor(temp[temp[,2]>1,1],temp[temp[,2]>1,2]))
    lm.fit(x = as.matrix(temp[temp[,2]>1,2]), temp[temp[,2]>1,1])$coefficients
  }, spp_mat, trait_dist)
  p_val <- list(lower=mean(temp[-1]<temp[1]), double = mean(abs(temp[-1]-mean(temp[-1]))>abs(temp[1]-mean(temp[-1]))))
  return(list(observed = temp[1], permutations = temp[-1], p = p_val))

} # TODO: This is not a good implementation, just works, you gotta change this later, in C++ or get it better in R


