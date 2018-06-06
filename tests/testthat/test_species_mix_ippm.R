context('species_mix ippm')

testthat::test_that('species mix functions classes work', {

  library(ecomix)
  library(raster)
  set.seed(42)

  x1 <- sort(runif(1000,0,2.5))
  x2 <- I(x1)^2

  n_g <- 4
  set.seed(123)
  thetas <- matrix(c(-7,38,-35,
                     -1.4,1.8,0,
                     -13,23,-8.2,
                     3.5,.2,-10.8),4,3,byrow=TRUE)

  dat <- data.frame(y=rep(1,100),x1,x2)

  x <- y <- 1:100 / 100
  grid2D <- expand.grid( x, y)
  grid2D$cellArea <- rep( 1/200, nrow( grid2D))  #all cells have same size here
  colnames(grid2D) <- c("x","y","cellArea")

  # now let's set up a variable to model.
  set.seed(6)
  d <- as.matrix(dist(grid2D[,c("x","y")]))
  w <- exp(-1/nrow(grid2D[,c("x","y")]) * d)
  ww <- chol(w)
  grid2D$x1 <- t(ww) %*% rnorm(nrow(grid2D[,c("x","y")]),0, 0.1)
  grid2D$x1 <- scales::rescale(grid2D$x1,to=range(0,2.5))

  coordinates(grid2D) <- ~x+y
  env <- rasterize(grid2D, raster(points2grid(grid2D)), fields=c("x1"))
  env<-dropLayer(env,1:2)

  n_sp <- 50
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  sp_name <- LETTERS702[1:(n_sp)]
  n_g <- 4

  # set.seed(123)
  X <- as.matrix(data.frame(const=1,x1=grid2D$x1,x2=I(grid2D$x1)^2))
  lambdas <- matrix(0, dim(X)[1], n_sp, dimnames=list(NULL,sp_name))
  sp_int <- rep(0, n_sp)
  group <- rep(0, n_sp)
  for (s in 1:n_sp) {
    g <- ceiling(runif(1) * n_g)
    sp_int[s] <- rnorm(1, thetas[g,1],.1)
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

  sp_weights_mat <- data.frame(cell_id = 1:10000, do.call(cbind,sp_weights))

  m <- sp_weights_mat
  presence_sites <- m[rowSums(is.na(m[,-1]))!=ncol(m[,-1]), ]
  presence_sites <- data.frame(presence_sites)

  background_sites <- data.frame(cell_id=1:10000,matrix(rep(grid2D$cellArea,n_sp),nrow(grid2D),n_sp))

  wts <- rbind(presence_sites[,-1],background_sites[,-1])
  colnames(wts) <- c(sp_name)
  wts <- as.matrix(wts)

  sam_form <- as.formula(paste0('cbind(',paste(LETTERS702[1:(n_sp)],collapse = ','),")~1+x1+x2"))
  sp_form <- ~ 1

  fm_ippm1 <- ecomix::species_mix(archetype_formula = sam_form, species_formula = sp_form, data = dat, weights = wts, distribution = 'ippm', n_mixtures = 4, titbits =  TRUE, control = species_mix.control(em_prefit = FALSE, calculate_hessian_cpp=FALSE))

  testthat::expect_s3_class(fm_ippm1,'species_mix')
  testthat::expect_s3_class(fm_ippm1,'ippm')
})
