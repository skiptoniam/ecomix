context('species_mix ippm')

testthat::test_that('species mix ippm', {

  library(ecomix)
  library(raster)
  set.seed(42)

  x1 <- sort(runif(1000,-2.5,2.5))
  x2 <- I(x1)^2

  n_g <- 4
  set.seed(123)
  thetas <- matrix(c( 1.6, 4.0, -3.0,
                      -0.8, 1.8,    0,
                      1.4,-1.2, -0.1,
                      1.2,-8.2, -6.8),4,3,byrow=TRUE)

  set.seed(42)
  x <- y <- 1:100 / 100
  grid2D <- expand.grid( x, y)
  grid2D$cellArea <- rep( 1/200, nrow( grid2D))  #all cells have same size here
  colnames(grid2D) <- c("x","y","cellArea")


  # now let's set up a variable to model.
  set.seed(6)
  d <- as.matrix(stats::dist(grid2D[,c("x","y")]))
  w <- exp(-1/nrow(grid2D[,c("x","y")]) * d)
  ww <- chol(w)
  grid2D$x1 <- t(ww) %*% rnorm(nrow(grid2D[,c("x","y")]),0, 0.1)
  grid2D$x1 <- scales::rescale(grid2D$x1,to=range(-2.5,2.5))

  coordinates(grid2D) <- ~x+y
  env <- rasterize(grid2D, raster(points2grid(grid2D)), fields=c("x1"))
  env<-dropLayer(env,1:2)

  n_sp <- 50
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  sp_name <- LETTERS702[1:(n_sp)]
  n_g <- 4

  set.seed(123)
  X <- as.matrix(data.frame(const=1,x1=grid2D$x1,x2=I(grid2D$x1)^2))
  lambdas <- matrix(0, dim(X)[1], n_sp, dimnames=list(NULL,sp_name))
  sp_int <- rep(0, n_sp)
  group <- rep(0, n_sp)
  for (s in 1:n_sp) {
    g <- sample(n_g,1)
    sp_int[s] <- rnorm(1, thetas[g,1], .5)
    log_lambda <-  X%*%c(sp_int[s],thetas[g,-1])
    lambdas[, s] <- exp(log_lambda)
    group[s] <- g
  }

  LAMBDAS <- apply(lambdas,2,function(x)sum(x*grid2D$cellArea))
  Ns <- sapply(LAMBDAS,function(x)rpois(n=1, lambda= x))  #the observed number of presences
  preds_df <- data.frame(idx=1:nrow(X),X)
  presences <- list()
  for(i in seq_len(n_sp)){
    presences[[i]] <- sample(x=preds_df$idx,size=Ns[i], replace=TRUE, prob=lambdas[,i]/LAMBDAS[i])# TRUE
  }

  presence_coords <- lapply(presences,function(x)coordinates(grid2D)[x,1:2])
  presences_sort <- lapply(presences,sort)

  sp_dat_po_ul<-data.frame(sp=rep(sp_name,unlist(lapply(presences,length))),cell_num=unlist(presences_sort))
  po_matrix <- bbgdm::table2pam(sp_dat_po_ul,site.id='cell_num',sp.id='sp')
  po_matrix[po_matrix==0]<-NA
  po_covariates <- X[as.numeric(rownames(po_matrix)),]
  presence_data <- data.frame(po_matrix,po_covariates)
  bkdata <- cbind(matrix(0,nrow(X),n_sp),X)
  colnames(bkdata) <- colnames(presence_data)
  mm <- rbind(presence_data,bkdata)
  dat <- mm[c(sp_name,"const","x1","x2")]

  species_specific_cell_counts <- lapply(seq_along(sp_name),function(x)table(sp_dat_po_ul[sp_dat_po_ul$sp==sp_name[x],2]))

  df <- data.frame(id=preds_df$idx,area=grid2D$cellArea,x1=grid2D$x1)

  sp_weights <- lapply(seq_along(sp_name),function(x)(weights=df$area/as.numeric(species_specific_cell_counts[[x]][match(df$id,as.numeric(names(species_specific_cell_counts[[x]])))])))

  # sp_weights_rows <- do.call(rbind,sp_weights)


  sp_weights_mat <- data.frame(cell_id = 1:10000, do.call(cbind,sp_weights))

  m <- sp_weights_mat
  presence_sites <- m[rowSums(is.na(m[,-1]))!=ncol(m[,-1]), ]
  presence_sites <- data.frame(presence_sites)#,const=1,x1=grid2D$x1[presence_sites$cell_id],x2=I(grid2D$x1[presence_sites$cell_id])^2)

  background_sites <- data.frame(cell_id=1:10000,matrix(rep(grid2D$cellArea,n_sp),nrow(grid2D),n_sp))#,const=1,x1=grid2D$x1,x2=I(grid2D$x1)^2)

  wts <- rbind(presence_sites[,-1],background_sites[,-1])
  colnames(wts) <- c(sp_name)#,"const","x1","x2")
  offset <- rep(0,nrow(dat))

  ## test the internal functions.
  ## test the apply functions for ippm

  y <- dat[,1:n_sp]
  X <- dat[,c(n_sp+1):ncol(dat)]
  y_is_na <- is.na(y)
  weights <- wts
  G <- 4

  # test if one species ippm working - expect matrix of coefs back
  one_sp_ippm <- ecomix:::apply_glmnet_ippm(1, y = y, X = X, weights = weights, offset = offset, y_is_na = y_is_na)
  testthat::expect_is(one_sp_ippm,'matrix')

  # check that many species ippms work - expect back a list.
  all_sp_ippm <-surveillance::plapply(1:n_sp, ecomix:::apply_glmnet_ippm, y, X, weights , offset , y_is_na )
  testthat::expect_is(all_sp_ippm,'list')

  #glmnet coefs
  all_coefs_mat <- t(do.call(cbind,all_sp_ippm))
  mix_coefs <- all_coefs_mat[,-1] # drop intercepts
  tmp1 <- kmeans(mix_coefs, centers=G, nstart=100)
  tmp_grp <- tmp1$cluster
  grp_coefs <- apply(mix_coefs, 2, function(x) tapply(x, tmp_grp, mean))

  # test if one species ippm working - expect matrix of coefs back
  one_sp_ippm_glm <- ecomix:::apply_glm_ippm(1, y = y, X = X, weights = weights, offset = offset, y_is_na = y_is_na)
  testthat::expect_is(one_sp_ippm,'matrix')

  # check that many species ippms work - expect back a list.
  all_sp_ippm_glm <-surveillance::plapply(1:n_sp, ecomix:::apply_glm_ippm, y, X, weights , offset , y_is_na )
  testthat::expect_is(all_sp_ippm,'list')

  ## make these into a matrix
  all_coefs_mat_glm <- t(do.call(cbind,all_sp_ippm_glm))
  mix_coefs_glm <- all_coefs_mat_glm[,-1] # drop intercepts

  tmp2 <- kmeans(mix_coefs_glm, centers=G, nstart=100)
  tmp_grp_glm <- tmp2$cluster
  grp_coefs_glm <- apply(mix_coefs_glm, 2, function(x) tapply(x, tmp_grp, mean))

  #now we need to estimate the taus.
  S <- 50
  G <- 4
  control <- species_mix.control()

  # expect error if wrong data is in the starting values
  testthat::expect_error(  starting_values <- ecomix:::initiate_fit_ippm(NULL, X, weights, offset, y_is_na, G, S, control))
  testthat::expect_error(  starting_values <- ecomix:::initiate_fit_ippm(y, NULL, weights, offset, y_is_na, G, S, control))
  testthat::expect_error(  starting_values <- ecomix:::initiate_fit_ippm(y, X, NULL, offset, y_is_na, G, S, control))
  testthat::expect_error(  starting_values <- ecomix:::initiate_fit_ippm(y, X, weights, offset, NULL, G, S, control))

  #expect list back
  starting_values <- ecomix:::initiate_fit_ippm(y, X, weights, offset, y_is_na, G, S, control)
  testthat::expect_is(starting_values,'list')

  fits <- list(betas=starting_values$mix_coefs, alphas=starting_values$sp_intercepts)
  first_fit <- list(x = X, y = y, weights=weights, offset=offset, y_is_na=y_is_na)

  # get the loglikelihood based on these values
  logls <- ecomix:::get_logls_ippm(first_fit, fits, G, S)
  testthat::expect_is(logls,'matrix')
  testthat::expect_equal(ncol(logls), G)
  testthat::expect_equal(nrow(logls), S)

  # estimate the posteriors for taus
  pis <- rep(1/G, G)
  taus <- ecomix:::get_taus_ippm(pis, logls, G, S)
  testthat::expect_is(taus,'matrix')
  testthat::expect_equal(ncol(taus), G)
  testthat::expect_equal(nrow(taus), S)

  # skrink the taus
  taus <- skrink_taus_ippm(taus, max_tau=0.99, G)
#

  ## now test if the group_tau glm works

  testthat::expect_error(fm_g1 <- ecomix:::apply_glmnet_ippm_group_tau_v2(gg = 1, y = y, X = X, weights = weights, offset = offset,
                                          y_is_na = y_is_na, tau = tau, return_all_coefs = FALSE))
  testthat::expect_error(fm_g1 <- ecomix:::apply_glmnet_ippm_group_tau_v2(gg = 1, y = NULL, X = X, weights = weights, offset = offset,
                                                   y_is_na = y_is_na, tau = taus, return_all_coefs = FALSE))
  testthat::expect_error(fm_g1 <- ecomix:::apply_glmnet_ippm_group_tau_v2(gg = 1, y = y, X = NULL, weights = weights, offset = 'blah',
                                                   y_is_na = y_is_na, tau = taus, return_all_coefs = FALSE))
  testthat::expect_error(fm_g1 <- ecomix:::apply_glmnet_ippm_group_tau_v2(gg = 1, y = y, X = X, weights = 'a', offset = offset,
                                                   y_is_na = y_is_na, tau = taus, return_all_coefs = FALSE))
  fm_g1 <- ecomix:::apply_glmnet_ippm_group_tau_v2(gg = 1, y = y, X = X, weights = weights, offset = offset,
                                                   y_is_na = y_is_na, tau = taus, return_all_coefs = FALSE)

  testthat::expect_is(fm_g1,'matrix')

  #now we test if it works with lapply
  all_grp_ippm1 <- surveillance::plapply(seq_len(G), ecomix:::apply_glmnet_ippm_group_tau_v2, y, X, weights, offset, y_is_na, taus)
  testthat::expect_is(all_grp_ippm1,'list')
  all_grp_ippm2 <- surveillance::plapply(seq_len(G), ecomix:::apply_glmnet_ippm_group_tau_v2, y, X, weights, offset, y_is_na, taus, lambda_pen = 1/2)
  testthat::expect_is(all_grp_ippm2,'list')
  all_grp_ippm3 <- surveillance::plapply(seq_len(G), ecomix:::apply_glmnet_ippm_group_tau_v2, y, X, weights, offset, y_is_na, taus, lambda_pen = 1/3)
  testthat::expect_is(all_grp_ippm3,'list')
  all_grp_ippm4 <- surveillance::plapply(seq_len(G), ecomix:::apply_glmnet_ippm_group_tau_v2, y, X, weights, offset, y_is_na, taus, lambda_pen = 1/5)
  testthat::expect_is(all_grp_ippm4,'list')
  all_grp_ippm5 <- surveillance::plapply(seq_len(G), ecomix:::apply_glmnet_ippm_group_tau_v2, y, X, weights, offset, y_is_na, taus, lambda_pen = 1/10)
  mix_coefs <- t(do.call(cbind,all_grp_ippm1))[,-1]


  # sam_form <- as.formula(paste0('cbind(',paste(LETTERS702[1:(n_sp)],collapse = ','),")~1+x1+x2"))
  # sp_form <- ~ 1
  #
  # fm_ippm1 <- ecomix::species_mix(archetype_formula = sam_form, species_formula = sp_form, data = dat, weights = as.matrix(weights),
  #                                 distribution = 'ippm', n_mixtures = 4, titbits =  TRUE,
  #                                 control = species_mix.control(em_prefit = FALSE, calculate_hessian_cpp=FALSE))
  #
  # testthat::expect_s3_class(fm_ippm1,'species_mix')
  # testthat::expect_s3_class(fm_ippm1,'ippm')
})
