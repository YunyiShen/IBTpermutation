make_4_plots <- function(res, filename, beta_size,
                         beta_dist, betas_env = c(0,-.5,-1,-3),...){

  pdf(filename,...)
  par(mfrow = c(2,2))
  labs <- c("A","B","C","D")
  for(i in 1:length( betas_env)){
    beta_env <- betas_env[i]
    p_vals <- res$p[res[,2]==beta_size & res[,3]==beta_dist & res[,4]==beta_env]
    plot(ecdf(p_vals), main = "", xlab = "p", ylab = "Fn(p)", xlim = c(0,1), ylim = c(0,1))
    txt <- labs[i]
    text(1 - strwidth(txt, cex=1.5) / 2, 0.05 + strheight(txt, cex=1.5) / 2, txt, cex=1.5)
    abline(0,1, col = "red")
    abline(v = .05, col = "blue")
  }

  dev.off()
}

p_vals_res <- read.csv("./tests/formal_simulations/logistic_IBT/p_vals_less_traits2.csv")

betas_size <- unique(p_vals_res$beta_size)
betas_dist <- unique(p_vals_res$beta_dist)
betas_env <- unique(p_vals_res$beta_env)

all_IBT_betas <- expand.grid(betas_size, betas_dist)

for(i in 1:nrow(all_IBT_betas) ){
  the_beta <- all_IBT_betas[i,]
  file_name <- paste0("./tests/formal_simulations/logistic_IBT/sorting_less_trait_", paste0(the_beta, collapse = "_"),".pdf")
  make_4_plots(p_vals_res, file_name, the_beta[1,1], the_beta[1,2],
               betas_env, width = 6, height = 6)
}


