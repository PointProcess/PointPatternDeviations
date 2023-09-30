# simETASspTemp <- function(mu, K, alpha, c, p, d, q, beta, M0, maxT, S.aux, S.target, displayOutput=TRUE, random.seed = NA, 
#                           mod.time=NA, ashock.duration = "fixed") {
#   # mu <- 0.1/100; K <- 0.006; alpha <- 2.2999; c <- 0.05; p <- 1.08; beta <- 2.3; d <- 1; q <- 2
#   # S.aux <- 1000; S.target <-500
#   # M0 <- 2.5; maxT <- 5000*100 ; random.seed <- NA; mod.time <- 5000; ashock.duration = "fixed"
#   # cat3.st <- simETASspTemp(mu=0.1/100, K=0.006, alpha=2.2999, c=0.05, p=1.08, d=1, q=2, beta=2.3, M0=2.5, 
#   # displayOutput = T, maxT=maxT*100, S.aux=1000, S.target=500, mod.time=maxT)
#   
#   # mu <- 0.1/10; K <- 0.9; alpha <- 1.5; c <- 0.5; p <- 1.04; beta <- 2.3
#   # M0 <- 2.5; maxT <- 2000*10 ;  mod.time <- 2000; ashock.duration = "varying"
#   # displayOutput = TRUE
#   
#   cp <-c
#   # Draw the number of background events ~ Pois(mu * length of time period)
#   if(!is.na(random.seed)){
#     set.seed(random.seed)
#   }
#   num.bkgd <- rpois(n=1, lambda=mu*maxT)
#   
#   # Times of bkgd events are uniformly distributed within the time period; 
#   # their magnitudes are GR distributed, using given beta
#   if(!is.na(random.seed)){
#     set.seed(random.seed)
#   }
#   # Add random spatial locations, uniformly distributed within 0 to S
#   bkgd.eqks <- data.frame("time" = runif(n=num.bkgd, min=0, max=maxT),
#                           "x" = runif(n=num.bkgd, min=0, max=S.aux),
#                           "y" = runif(n=num.bkgd, min=0, max=S.aux),
#                           "mag"= GRFunc(n=num.bkgd, beta=beta, M0=M0))
#   # TODO: Remove this. I mod T all the times later on
#   # if(!is.na(mod.time)){
#   # bkgd.eqks$time <- bkgd.eqks$time %% mod.time
#   # }
#   bkgd.eqks <- bkgd.eqks[order(bkgd.eqks$time),]
#   bkgd.eqks$index <- 1:num.bkgd
#   bkgd.eqks$gen <- 0
#   bkgd.eqks$branching <- 0
#   all.eqks <- bkgd.eqks
#   
#   i <- 1
#   # While there are still eqks that can trigger others...
#   # Take eqk i...
#   while (i < nrow(all.eqks)){
#     # The aftershock period duration should be 1 or 2x the unscaled catalog period
#     if(ashock.duration == "fixed"){
#       maxT.trig <- maxT*2  # Can also try 2x the catalog period
#       # So if we're using mod T procedure, maxT.trig is mod.time
#       if(!is.na(mod.time)){
#         maxT.trig <- mod.time*2 
#       }
#     } else if(ashock.duration == "varying"){
#       maxT.trig <- maxT-all.eqks$time[i] 
#       # So if we're using mod T procedure, maxT.trig is remaining time from desired cat length (mod.time)
#       if(!is.na(mod.time)){
#         maxT.trig <-  mod.time-(all.eqks$time[i] %% mod.time)
#       }
#     }
#     # remaining.T <- maxT - all.eqks$time[i]
#     
#     # Integrate the Omori function over the time interval (0, maxT.trig) and scale by the productivity 
#     # to get the mean number of events expected to be triggered, by this eqk. The number triggered is 
#     # based SOLELY on the mainshock's mag under a "fixed aftershock duration" strategy
#     mean.num.trig <- K*exp(alpha*(all.eqks$mag[i]-M0))*
#       ((cp+maxT.trig)^(1-p) - cp^(1-p))*(1/(1-p)) # Before, the last term was -cp^(p-1)
#     # mean.num.trig <- K*exp(alpha*(all.eqks$mag[i]-M0))*
#     #   ((cp+remaining.T)^(1-p) - cp^(1-p))/(1-p)
#     # print(all.eqks$time[i])
#     all.eqks$num.trig[i] <- rpois(n=1, lambda=mean.num.trig)
#     
#     print(i) ; print(all.eqks$num.trig[i])
#     # If any events were triggered...
#     if(all.eqks$num.trig[i] > 0){
#       # Draw the times of these triggered events according to the Omori distribution
#       # Draw the magnitudes of these triggered events according to the GR distribution
#       
#       # print(paste0("Cat length is ", maxT, "and aftershock max length is ", maxT.trig))
#       trig.eqks <- newGenTrig(current.eqk=all.eqks[i,], beta=beta, current.eqk.i=i, 
#                               maxT.trig=maxT.trig, maxT.cat=maxT, M0=M0, cp=cp, p=p, d=d, q=q)
#       # random.seed = random.seed, 
#       # Add index values starting from the last eqk in the set and bind these triggered eqks
#       # to the df of all eqks
#       if(nrow(trig.eqks) > 0){
#         trig.eqks$index <- nrow(all.eqks) + 1:nrow(trig.eqks)
#       }
#       # if(all(trig.eqks$gen < 3)){
#       all.eqks <- rbind(all.eqks, trig.eqks)
#       # }
#     }
#     # Move on to the next eqk in the df
#     i <- i + 1
#   }
#   
#   # Sort chronologically
#   all.eqks <- all.eqks[order(all.eqks$time),]
#   # Calculate time diffs
#   all.eqks$time.diff <- calcTimeDiffs(all.eqks)
#   
#   # If we're mod-timing..
#   if(!is.na(mod.time)){
#     print("Using transformed times!")
#     
#     # TODO: Remove? 
#     # Sort by gen to avoid fallen seqs within fallen seqs
#     all.eqks <- all.eqks[order(all.eqks$gen),]
#     all.eqks$mod.time <- all.eqks$time %% mod.time
#     
#     # Find all triggered events
#     trig.events.all <- all.eqks$index[all.eqks$branching !=0]
#     eqks.falling.out.of.mod.time <- c()
#     # For each triggered event...
#     for(i in trig.events.all){
#       # Pull out its time and the time of the event that triggered it 
#       # i <- trig.events.all[5]
#       this.time <- all.eqks$mod.time[all.eqks$index == i]
#       this.trigger <- all.eqks$branching[all.eqks$index == i]
#       this.triggers.time <- all.eqks$mod.time[all.eqks$index == this.trigger]
#       # print(this.triggers.time - this.time)
#       # If its trigger's mod time is greater than its mod time, its fallen out of its seq 
#       if(this.triggers.time > this.time){
#         eqks.falling.out.of.mod.time <- c(i, eqks.falling.out.of.mod.time)
#       }
#     }
#     print(length(eqks.falling.out.of.mod.time))
#     
#     # all.eqks.og <- all.eqks
#     # all.eqks <- all.eqks.og
#     
#     # Find the new branching structure of the catalog cut by mod time
#     all.eqks$cut.by.mod <- FALSE
#     all.eqks$mod.gen <- all.eqks$gen 
#     all.eqks$mod.branching <- all.eqks$branching 
#     # For each new bkgd event (that separated from its direct trigger)
#     for(this.eqk in eqks.falling.out.of.mod.time){
#       # print(paste0("this eqk is ", this.eqk))
#       # this.eqk <-  488
#       # Find its gen and branching tree. All of these were cut by the mod procedure
#       gen.this.eqk <- all.eqks$gen[all.eqks$index == this.eqk]
#       
#       # all.eqks[all.eqks$index == this.eqk,]
#       # all.eqks[all.eqks$index == 394,]
#       
#       this.eqk.full.tree <- eqkTreeAllBranches(all.eqks, this.eqk)
#       # all.eqks[all.eqks$index %in% this.eqk.full.tree,]
#       # this.eqk.full.tree.dt <- all.eqks[all.eqks$index %in% this.eqk.full.tree,]
#       
#       # all.eqks[all.eqks$index %in% this.eqk.full.tree,]
#       all.eqks$cut.by.mod[all.eqks$index %in% this.eqk.full.tree] <- TRUE
#       
#       for(this.trig.eqk in unique(this.eqk.full.tree)){
#         # print(this.trig.eqk)
#         # If we haven't already changed this eqk's gen...
#         # print(all.eqks$gen[all.eqks$index == this.trig.eqk])
#         if(all.eqks$mod.gen[all.eqks$index == this.trig.eqk] == 
#            all.eqks$gen[all.eqks$index == this.trig.eqk]){
#           # Then reduce it by the gen of the eqk that was dropped to gen 0
#           all.eqks$mod.gen[all.eqks$index == this.trig.eqk] <-
#             all.eqks$gen[all.eqks$index == this.trig.eqk] - gen.this.eqk
#         }
#         # print(all.eqks$mod.gen[all.eqks$index == this.trig.eqk])
#       }
#       
#       # Eqks that were split by mod time from their mainshock should get branching=0 
#       # (they are now new bkgd events)
#       all.eqks$mod.branching[all.eqks$index == this.eqk] <- 0
#       
#     } 
#     
#     # Check that no gens are negative
#     stopifnot(all(all.eqks$mod.gen >= 0))
#     
#     # All those whose new gen is 0 should be either original bkgd events or those
#     # that were cut off
#     stopifnot(all(all.eqks$index[(all.eqks$new.gen == 0)] %in% 
#                     c(all.eqks$index[all.eqks$branching == 0], 
#                       eqks.falling.out.of.mod.time)))
#     
#     stopifnot(all(all.eqks$index[(all.eqks$new.branching == 0)] %in% 
#                     c(all.eqks$index[all.eqks$branching == 0], 
#                       eqks.falling.out.of.mod.time)))
#     
#     stopifnot(all(all.eqks$index[all.eqks$new.gen != all.eqks$gen] %in% 
#                     all.eqks$index[all.eqks$cut.by.mod]))
#     
#     # bkgdPropCat(all.eqks)
#     # Change time to be the transformed time and call original time raw.time
#     all.eqks$time.raw <- all.eqks$time
#     all.eqks$time <- all.eqks$mod.time
#     # Remove the mod.time var
#     all.eqks$mod.time <- NULL
#     # Sort chronologically by mod time
#     all.eqks <- all.eqks[order(all.eqks$time),]
#     # Calculate time diffs using the new branching structure
#     all.eqks$time.diff <- calcTimeDiffs(all.eqks, branching.var = "mod.branching")
#   }
#   
#   # Spatial boundary condition: only include points within the target window
#   midpoint <- (S.aux - S.target)/2
#   all.eqks$target <- FALSE
#   all.eqks$target[all.eqks$x > midpoint & all.eqks$x < S.target + midpoint &
#                     all.eqks$y > midpoint & all.eqks$y < S.target + midpoint] <- TRUE
#   # plot(all.eqks$x[all.eqks$target], all.eqks$y[all.eqks$target])  
#   # table(all.eqks$target)
#   
#   return(all.eqks)
# }

runSTBETAS <- function(cat, init.vals, Geval.mat, fixed.br=F, S.target=NA, spatial.area=NA, sims=5000, maxT=20000, isAlphaFixed=F, M0){
  # cat <- pnw.Saux.comp1.cat.noCswarms
  # init.vals <- cat2inits
  # Geval.mat <- pnw.Saux.comp1.cat.noCswarms.Gevals.mat
  # fixed.br=F; S.target=NA; spatial.area=intl.Saux.area
  # sims=5000;  maxT=maxT.pnw.Saux; 
  # isAlphaFixed=F; M0 <- 2.3
  
  mu.cat <- init.vals[1]; K.cat <- init.vals[2]; alpha.cat <- init.vals[3]; 
  cp.cat <- init.vals[4]; p.cat <- init.vals[5]; d.cat <- init.vals[6]; q.cat <- init.vals[7];
  
  if(!fixed.br){
    print("YES!!")
    # out.cat.target.G.interp <- STestimateETASBranchingInteractionMSList(ts=as.double(cat$time),
    #                                                                     marks=as.double(cat$mag), lons=as.double(cat$x),
    #                                                                     lats=as.double(cat$y), 
    #                                                                     # branching=as.integer(cat$reindex.branching),
    #                                                                     branching=rep(0, length(cat$time)),
    #                                                                     maxT = as.double(maxT), 
    #                                                                     M0=as.double(2.5), sims=as.integer(sims), numMCMCSamples=as.integer(numMCMCSamples),
    #                                                                     mu=as.double(mu.cat), logK=as.double(log(K.cat)),
    #                                                                     alpha=as.double(alpha.cat), c=as.double(cp.cat), p=as.double(p.cat),
    #                                                                     d=as.double(d.cat), q=as.double(q.cat),
    #                                                                     cat_Gevals_grid = Geval.mat,
    #                                                                     mus=mus, logKs=logKs, alphas=alphas,
    #                                                                     cs=cs, ps=ps,
    #                                                                     ds=ds, qs=qs)
    if(!is.na(S.target) & is.na(spatial.area)){
      out.cat.target.G.interp <- STestimateETASBranchingFreeBMSList(ts=as.double(cat$time),
                                                                    marks=as.double(cat$mag), lons=as.double(cat$x),
                                                                    lats=as.double(cat$y), 
                                                                    # branching=as.integer(cat$reindex.branching),
                                                                    branching=rep(0, length(cat$time)),
                                                                    maxT = as.double(maxT), 
                                                                    M0=as.double(M0), sims=as.integer(sims), numMCMCSamples=as.integer(numMCMCSamples),
                                                                    mu=as.double(mu.cat), logK=as.double(log(K.cat)),
                                                                    alpha=as.double(alpha.cat), c=as.double(cp.cat), p=as.double(p.cat),
                                                                    d=as.double(d.cat), q=as.double(q.cat),
                                                                    initval = init.vals,
                                                                    cat_Gevals_grid = Geval.mat,
                                                                    Starget = S.target, alphaFixed = as.logical(isAlphaFixed))
      
    } else if(is.na(S.target) & !is.na(spatial.area)){
      out.cat.target.G.interp <- STestimateETASBranchingFreeBMSList(ts=as.double(cat$time),
                                                                    marks=as.double(cat$mag), lons=as.double(cat$x),
                                                                    lats=as.double(cat$y), 
                                                                    # branching=as.integer(cat$reindex.branching),
                                                                    branching=rep(0, length(cat$time)),
                                                                    maxT = as.double(maxT), 
                                                                    M0=as.double(M0), sims=as.integer(sims), numMCMCSamples=as.integer(numMCMCSamples),
                                                                    mu=as.double(mu.cat), logK=as.double(log(K.cat)),
                                                                    alpha=as.double(alpha.cat), c=as.double(cp.cat), p=as.double(p.cat),
                                                                    d=as.double(d.cat), q=as.double(q.cat),
                                                                    initval = init.vals,
                                                                    cat_Gevals_grid = Geval.mat,
                                                                    spatialArea = spatial.area, alphaFixed = as.logical(isAlphaFixed))
      
    }
    
  } else{
    # out.cat.target.G.interp <- STestimateETASBranchingFixBMSList(ts=as.double(cat$time),
    #                                                                     marks=as.double(cat$mag), lons=as.double(cat$x),
    #                                                                     lats=as.double(cat$y), 
    #                                                                     branching=as.integer(cat$reindex.branching),
    #                                                                     # branching=rep(0, length(cat$time)),
    #                                                                     maxT = as.double(maxT), 
    #                                                                     M0=as.double(2.5), sims=as.integer(sims), numMCMCSamples=as.integer(numMCMCSamples),
    #                                                                     mu=as.double(mu.cat), logK=as.double(log(K.cat)),
    #                                                                     alpha=as.double(alpha.cat), c=as.double(cp.cat), p=as.double(p.cat),
    #                                                                     d=as.double(d.cat), q=as.double(q.cat),
    #                                                                     cat_Gevals_grid = Geval.mat,
    #                                                                     mus=mus, logKs=logKs, alphas=alphas,
    #                                                                     cs=cs, ps=ps,
    #                                                                     ds=ds, qs=qs)
    out.cat.target.G.interp <- STestimateETASBranchingFixBMSList(ts=as.double(cat$time),
                                                                 marks=as.double(cat$mag), lons=as.double(cat$x),
                                                                 lats=as.double(cat$y), 
                                                                 branching=as.integer(cat$reindex.branching),
                                                                 # branching=rep(0, length(cat$time)),
                                                                 maxT = as.double(maxT), 
                                                                 M0=as.double(M0), sims=as.integer(sims), numMCMCSamples=as.integer(numMCMCSamples),
                                                                 mu=as.double(mu.cat), logK=as.double(log(K.cat)),
                                                                 alpha=as.double(alpha.cat), c=as.double(cp.cat), p=as.double(p.cat),
                                                                 d=as.double(d.cat), q=as.double(q.cat),
                                                                 initval = init.vals,
                                                                 cat_Gevals_grid = Geval.mat)
  }
  posts.cat.target.G.interp <- as.data.frame(cbind(out.cat.target.G.interp$nbkgd, out.cat.target.G.interp$mu, 
                                                   out.cat.target.G.interp$K, out.cat.target.G.interp$alpha,
                                                   out.cat.target.G.interp$c, out.cat.target.G.interp$p, out.cat.target.G.interp$d, 
                                                   out.cat.target.G.interp$q))
  # posts.cat.target.G.interp  <- posts.cat.target.G.interp[-c(1:burnin, nrow(posts.cat.target.G.interp)),]
  
  names(posts.cat.target.G.interp) <- c("nbkgd","mu", "K", "alpha", "c", "p", "d", "q")
  # posts.cat.target.G.interp$mu <- posts.cat.target.G.interp$mu*(500^2)
  # posts.cat.target.G.interp$mu <- posts.cat.target.G.interp$mu*500/((500^2)/(1000^2))
  
  return(posts.cat.target.G.interp)
  
}


runSTBETASPrior <- function(cat, init.vals, Geval.mat, S.target=NA, spatial.area=NA, 
                            sims=5000, maxT=20000, isAlphaFixed=F, M0, pgammas_vec,
                            mu.hypers, K.hypers, alpha.hypers, this.M0, 
                            c.hypers, p.hypers, d.hypers, q.hypers){
  mu.cat <- init.vals[1]; K.cat <- init.vals[2]; alpha.cat <- init.vals[3]; 
  cp.cat <- init.vals[4]; p.cat <- init.vals[5]; d.cat <- init.vals[6]; q.cat <- init.vals[7];
  
  if(!is.na(S.target) & is.na(spatial.area)){
    out.cat.target.G.interp <- STestimateETASBranchingFreeBMSListPrior(ts=as.double(cat$time),
                                                                       marks=as.double(cat$mag), lons=as.double(cat$x),
                                                                       lats=as.double(cat$y), 
                                                                       # branching=as.integer(cat$reindex.branching),
                                                                       branching=rep(0, length(cat$time)),
                                                                       maxT = as.double(maxT), 
                                                                       M0=as.double(this.M0), sims=as.integer(sims), numMCMCSamples=as.integer(numMCMCSamples),
                                                                       mu=as.double(mu.cat), logK=as.double(log(K.cat)),
                                                                       alpha=as.double(alpha.cat), c=as.double(cp.cat), p=as.double(p.cat),
                                                                       d=as.double(d.cat), q=as.double(q.cat),
                                                                       initval = init.vals,
                                                                       cat_Gevals_grid = Geval.mat,
                                                                       Starget = S.target, 
                                                                       alphaFixed = as.logical(isAlphaFixed),
                                                                       mupgamma = as.logical(pgammas_vec[1]), Kpgamma = as.logical(pgammas_vec[2]),
                                                                       alphapgamma = as.logical(pgammas_vec[3]), cpgamma = as.logical(pgammas_vec[4]),
                                                                       ppgamma = as.logical(pgammas_vec[5]), dpgamma = as.logical(pgammas_vec[6]),
                                                                       qpgamma = as.logical(pgammas_vec[7]),
                                                                       muhyper1 = mu.hypers[1], muhyper2 = mu.hypers[2],
                                                                       Khyper1 = K.hypers[1], Khyper2 = K.hypers[2],
                                                                       alphahyper1 = alpha.hypers[1], alphahyper2 = alpha.hypers[2],
                                                                       chyper1 = c.hypers[1], chyper2 = c.hypers[2],
                                                                       phyper1 = p.hypers[1], phyper2 = p.hypers[2],
                                                                       dhyper1 = d.hypers[1], dhyper2 = d.hypers[2],
                                                                       qhyper1 = q.hypers[1], qhyper2 = q.hypers[2])
    
  } else if(is.na(S.target) & !is.na(spatial.area)){
    out.cat.target.G.interp <- STestimateETASBranchingFreeBMSListPrior(ts=as.double(cat$time),
                                                                       marks=as.double(cat$mag), lons=as.double(cat$x),
                                                                       lats=as.double(cat$y), 
                                                                       # branching=as.integer(cat$reindex.branching),
                                                                       branching=rep(0, length(cat$time)),
                                                                       maxT = as.double(maxT), 
                                                                       M0=as.double(this.M0), sims=as.integer(sims), numMCMCSamples=as.integer(numMCMCSamples),
                                                                       mu=as.double(mu.cat), logK=as.double(log(K.cat)),
                                                                       alpha=as.double(alpha.cat), c=as.double(cp.cat), p=as.double(p.cat),
                                                                       d=as.double(d.cat), q=as.double(q.cat),
                                                                       initval = init.vals,
                                                                       cat_Gevals_grid = Geval.mat,
                                                                       spatialArea = spatial.area, alphaFixed = as.logical(isAlphaFixed),
                                                                       mupgamma = as.logical(pgammas_vec[1]), Kpgamma = as.logical(pgammas_vec[2]),
                                                                       alphapgamma = as.logical(pgammas_vec[3]), cpgamma = as.logical(pgammas_vec[4]),
                                                                       ppgamma = as.logical(pgammas_vec[5]), dpgamma = as.logical(pgammas_vec[6]),
                                                                       qpgamma = as.logical(pgammas_vec[7]),
                                                                       muhyper1 = mu.hypers[1], muhyper2 = mu.hypers[2],
                                                                       Khyper1 = K.hypers[1], Khyper2 = K.hypers[2],
                                                                       alphahyper1 = alpha.hypers[1], alphahyper2 = alpha.hypers[2],
                                                                       chyper1 = c.hypers[1], chyper2 = c.hypers[2],
                                                                       phyper1 = p.hypers[1], phyper2 = p.hypers[2],
                                                                       dhyper1 = d.hypers[1], dhyper2 = d.hypers[2],
                                                                       qhyper1 = q.hypers[1], qhyper2 = q.hypers[2])
    
  }
  
  posts.cat.target.G.interp <- as.data.frame(cbind(out.cat.target.G.interp$nbkgd, out.cat.target.G.interp$mu, 
                                                   out.cat.target.G.interp$K, out.cat.target.G.interp$alpha,
                                                   out.cat.target.G.interp$c, out.cat.target.G.interp$p, out.cat.target.G.interp$d, 
                                                   out.cat.target.G.interp$q))
  # posts.cat.target.G.interp  <- posts.cat.target.G.interp[-c(1:burnin, nrow(posts.cat.target.G.interp)),]
  
  names(posts.cat.target.G.interp) <- c("nbkgd","mu", "K", "alpha", "c", "p", "d", "q")
  # posts.cat.target.G.interp$mu <- posts.cat.target.G.interp$mu*(500^2)
  # posts.cat.target.G.interp$mu <- posts.cat.target.G.interp$mu*500/((500^2)/(1000^2))
  
  return(posts.cat.target.G.interp)
  
}

# runSTBETASPlots(cat=cat3A.nomod[cat3A.nomod$target & cat3A.nomod$time > 5000,], 
#                 isCat2=F, isTarget=T, cat.label="cat3A-nomod-target-tempbdry-initstrue", 
#                 cat.Gevals.mat=cat3A.Gevals.mat, this.S.target=c(250,750), 
#                 true.vals=cat3A.true, these.init.vals = abs(perturbInits(cat3A.true)), 
#                 plot.title = "Cat3A No mod (Temp Bdry, T>5000), Target, Init Vals Near True", 
#                 nsims=5000+burnin, maxT=15000)
runSTBETASPlots <- function(cat, isCat2, isTarget, cat.label, 
                            cat.Gevals.mat, this.S.target=c(250,750), true.vals, 
                            these.init.vals, plot.title, nsims, maxT=20000) {
  # cat=cat2D.nomod[cat2D.nomod$target,]; isCat2=T; isTarget=T; cat.label="cat2D-test"
  # cat.Gevals.mat=cat2D.Gevals.mat; this.S.target=c(250,750)
  # true.vals=cat2D.true; these.init.vals = abs(perturbInits(cat2D.true))
  # plot.title = "Cat2D Test"
  # nsims=5+burnin; maxT=15000
  # 
  # cat=intl.t; isCat2=F; isTarget=T; cat.label="pnw"
  # cat.Gevals.mat=cat2D.Gevals; this.S.target=c(250,750)
  # true.vals=cat2D.true; these.init.vals = abs(perturbInits(cat2D.true))
  # plot.title = "Cat3B No mod, Target, Init Vals Near True"
  # nsims=5+burnin; maxT=15000
  # print("weeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee")
  
  cat.posts <- runSTBETAS(cat, init.vals=these.init.vals, 
                          Geval.mat=cat.Gevals.mat, 
                          fixed.br=F, S.target=this.S.target, sims=nsims, maxT=maxT)
  if(isCat2){
    cat.string <- "cat2"
  } else{
    cat.string <- "cat3"
  }
  write.csv(cat.posts, 
            paste0("data/bayesianETAS-output/2020-experiments/st-model/", cat.string,
                   "/posts-", cat.label, ".csv"), row.names=F)
  cat.posts.plot <- cat.posts[-c(1:burnin),]
  if(isTarget){
    cat.posts.plot$mu <- cat.posts.plot$mu * (1000^2) * 4
    all.true.vals = c(sum(cat$reindex.branching == 0 & cat$target), true.vals)
  } else {
    cat.posts.plot$mu <- cat.posts.plot$mu * (1000^2)
    all.true.vals = c(sum(cat$reindex.branching == 0), true.vals)
  }
  
  
  plotPosteriorST(cat, cat.posts.plot, this.n=nsims, model = "spat-temp",
                  M0 = this.M0,
                  these.true.vals = all.true.vals,
                  plot.path= paste0("plots/models/bayesianETAS/2020-experiments/st-model/", cat.string, "/st-model-outputs/posts-", cat.label, ".pdf"))
  # plot.path= paste0("plots/models/bayesianETAS/2020-experiments/st-model/", cat.string, "/st-model-outputs/multi-realns/", substr(cat.label, 1, 5), "/posts-", cat.label, ".pdf"))
  plotTraceST(cat.posts.plot, 
              plot.path = paste0("plots/models/bayesianETAS/2020-experiments/st-model/", cat.string, "/st-model-outputs/trace-", cat.label, ".pdf"),
              # plot.path= paste0("plots/models/bayesianETAS/2020-experiments/st-model/", cat.string, "/st-model-outputs/multi-realns/", substr(cat.label, 1, 5), "/trace-", cat.label, ".pdf"),
              these.true.vals= all.true.vals)
  plotScatterMatrix(cat.posts.plot, plot.title = plot.title,
                    plot.path = paste0("plots/models/bayesianETAS/2020-experiments/st-model/", cat.string, "/st-model-outputs/smat-", cat.label, ".pdf"))
  # plot.path= paste0("plots/models/bayesianETAS/2020-experiments/st-model/", cat.string, "/st-model-outputs/multi-realns/", substr(cat.label, 1, 5), "/smat-", cat.label, ".pdf"))
  
}

runSTBETASPlotsPNW <- function(cat, cat.label, nsims, these.init.vals, 
                               this.spat.area = intl.t.area, this.isAlphaFixed=T,
                               cat.Gevals.mat, this.maxT, this.plot.loc,
                               plot.title, this.M0) {
  # cat=cat2D.nomod[cat2D.nomod$target,]; isCat2=T; isTarget=T; cat.label="cat2D-test"
  # cat.Gevals.mat=cat2D.Gevals.mat; this.S.target=c(250,750)
  # true.vals=cat2D.true; these.init.vals = abs(perturbInits(cat2D.true))
  # plot.title = "Cat2D Test"
  # nsims=5+burnin; maxT=15000
  # 
  # cat=intl.t; isCat2=F; isTarget=T; cat.label="pnw"
  # cat.Gevals.mat=cat2D.Gevals; this.S.target=c(250,750)
  # true.vals=cat2D.true; these.init.vals = abs(perturbInits(cat2D.true))
  # plot.title = "Cat3B No mod, Target, Init Vals Near True"
  # nsims=5+burnin; maxT=15000
  # print("weeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee")
  
  cat.posts <- runSTBETAS(cat, init.vals=these.init.vals, 
                          Geval.mat=cat.Gevals.mat, 
                          fixed.br=F, spatial.area=this.spat.area, sims=nsims, 
                          maxT=this.maxT, isAlphaFixed = this.isAlphaFixed, M0=this.M0)
  write.csv(cat.posts, 
            paste0("data/bayesianETAS-output/2020-experiments/st-model/pnw/posts-", cat.label, ".csv"), row.names=F)
  cat.posts.plot <- cat.posts[-c(1:burnin),]
  cat.posts.plot$mu <- cat.posts.plot$mu * this.spat.area
  plotPosteriorST(cat, cat.posts.plot, this.n=nsims, model = "spat-temp",
                  M0 = this.M0,
                  plot.path= paste0("plots/models/bayesianETAS/2020-experiments/st-model/", 
                                    this.plot.loc, "/posts-", cat.label, ".pdf"))
  # plot.path= paste0("plots/models/bayesianETAS/2020-experiments/st-model/", cat.string, "/multi-realns/", substr(cat.label, 1, 5), "/posts-", cat.label, ".pdf"))
  plotTraceST(cat.posts.plot, 
              plot.path = paste0("plots/models/bayesianETAS/2020-experiments/st-model/", 
                                 this.plot.loc, "/trace-", cat.label, ".pdf"))
  plotScatterMatrix(cat.posts.plot, plot.title = plot.title,
                    plot.path = paste0("plots/models/bayesianETAS/2020-experiments/st-model/", 
                                       this.plot.loc, "/smat-", cat.label, ".pdf"))
  # plot.path= paste0("plots/models/bayesianETAS/2020-experiments/st-model/", cat.string, "/multi-realns/", substr(cat.label, 1, 5), "/smat-", cat.label, ".pdf"))
  
}

runSTBETASBr <- function(cat, cat.label, nsims, these.init.vals, 
                               this.spat.area, this.isAlphaFixed=T,
                         Geval.mat, this.maxT, this.M0) {
  # cat=cat2D.nomod[cat2D.nomod$target,]; isCat2=T; isTarget=T; cat.label="cat2D-test"
  # cat.Gevals.mat=cat2D.Gevals.mat; this.S.target=c(250,750)
  # true.vals=cat2D.true; these.init.vals = abs(perturbInits(cat2D.true))
  # plot.title = "Cat2D Test"
  # nsims=5+burnin; maxT=15000
  # 
  # cat=pnw.Naux.cat.noCswarms; cat.label="pnw"
  # Geval.mat=pnw.Naux.cat.noCswarms.Gevals.mat; this.spat.area = intl.Naux.area
  # cat=intl.aux; cat.label="pnw"
  # Geval.mat=pnw.Gevals.aux.mat; this.spat.area = intl.Naux.area
  # these.init.vals = postmodealphafix.inits
  # plot.title = "Cat3B No mod, Target, Init Vals Near True"
  # nsims=1000; this.isAlphaFixed=T
  # this.maxT <- maxT.pnw.Naux; this.M0 <- 2.0
  print("weeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee")
  
  mu.cat <- these.init.vals[1]; K.cat <- these.init.vals[2]; alpha.cat <- these.init.vals[3]; 
  cp.cat <- these.init.vals[4]; p.cat <- these.init.vals[5]; d.cat <- these.init.vals[6]; q.cat <- these.init.vals[7];
  
  out.cat.target.G.interp <- STestimateETASBranchingFreeBMSList(ts=as.double(cat$time),
                                                                marks=as.double(cat$mag), lons=as.double(cat$x),
                                                                lats=as.double(cat$y), 
                                                                # branching=as.integer(cat$reindex.branching),
                                                                branching=rep(0, length(cat$time)),
                                                                maxT = as.double(this.maxT), 
                                                                M0=as.double(this.M0), sims=as.integer(nsims), numMCMCSamples=as.integer(numMCMCSamples),
                                                                mu=as.double(mu.cat), logK=as.double(log(K.cat)),
                                                                alpha=as.double(alpha.cat), c=as.double(cp.cat), p=as.double(p.cat),
                                                                d=as.double(d.cat), q=as.double(q.cat),
                                                                initval = these.init.vals,
                                                                cat_Gevals_grid = Geval.mat,
                                                                spatialArea = this.spat.area, alphaFixed = as.logical(this.isAlphaFixed))
  
  
  # posts.cat.target.G.interp <- as.data.frame(cbind(out.cat.target.G.interp$nbkgd, out.cat.target.G.interp$mu, 
  #                                                  out.cat.target.G.interp$K, out.cat.target.G.interp$alpha,
  #                                                  out.cat.target.G.interp$c, out.cat.target.G.interp$p, out.cat.target.G.interp$d, 
  #                                                  out.cat.target.G.interp$q))
  # # posts.cat.target.G.interp  <- posts.cat.target.G.interp[-c(1:burnin, nrow(posts.cat.target.G.interp)),]
  # 
  # names(posts.cat.target.G.interp) <- c("nbkgd","mu", "K", "alpha", "c", "p", "d", "q")
  # posts.cat.target.G.interp$mu <- posts.cat.target.G.interp$mu*(500^2)
  # posts.cat.target.G.interp$mu <- posts.cat.target.G.interp$mu*500/((500^2)/(1000^2))
  
  posts.br.cat.target.G.interp <- out.cat.target.G.interp$branching
  write.csv(posts.br.cat.target.G.interp, 
            paste0("data/bayesianETAS-output/2020-experiments/st-model/pnw/post-br-", cat.label, ".csv"), row.names=F)
  
  return(posts.br.cat.target.G.interp)
}

runSTBETASPlotsPNWPriors <- function(cat, cat.label, nsims, these.init.vals, this.M0,
                                     this.spat.area = intl.t.area, this.isAlphaFixed=T,
                                     cat.Gevals.mat, this.maxT, this.plot.loc, this.pgammas_vec,
                                     plot.title, these.mu.hypers, these.K.hypers, these.alpha.hypers, 
                                     these.c.hypers, these.p.hypers, these.d.hypers, these.q.hypers,
                                     this.PB.region = NA, this.JG.region = NA) {
  # cat=intl.t.noCswarms
  # cat.label = "pnw-noCswarms-Nt-cat2inits-fixalpha-JGpriors"
  # these.init.vals=cat2inits; this.M0=2.0
  # cat.Gevals.mat=pnw.t.noC.swarms.Gevals.mat
  # this.plot.loc="pnw/priors/"
  # this.spat.area=intl.t.area; nsims=5000
  # 
  # this.maxT=maxT.pnw; this.isAlphaFixed=T
  # these.mu.hypers=mu.SZ.Unif.hypers
  # these.K.hypers=K.SZ.Unif.hypers
  # these.alpha.hypers=alpha.SZ.Unif.hypers
  # these.c.hypers=c.SZ.Unif.hypers
  # these.p.hypers=p.SZ.Unif.hypers
  # these.d.hypers=d.SZ.Unif.hypers
  # these.q.hypers=q.SZ.Unif.hypers
  # this.pgammas_vec <- rep(FALSE, 7)
  
  cat.posts <- runSTBETASPrior(cat, init.vals=these.init.vals, 
                               Geval.mat=cat.Gevals.mat, this.M0=this.M0,
                               spatial.area=this.spat.area, sims=nsims, 
                               maxT=this.maxT, isAlphaFixed = this.isAlphaFixed, pgammas_vec=this.pgammas_vec,
                               mu.hypers=these.mu.hypers,
                               K.hypers=these.K.hypers,alpha.hypers=these.alpha.hypers,
                               c.hypers=these.c.hypers,p.hypers=these.p.hypers,
                               d.hypers=these.d.hypers,q.hypers=these.q.hypers)
  write.csv(cat.posts, 
            paste0("data/bayesianETAS-output/2020-experiments/st-model/pnw/posts-", cat.label, ".csv"), row.names=F)
  cat.posts.plot <- cat.posts[-c(1:burnin),]
  cat.posts.plot$mu <- cat.posts.plot$mu * this.spat.area
  
  this.hypers.list <- list(these.mu.hypers, these.K.hypers, these.alpha.hypers,
                     these.c.hypers, these.p.hypers,
                      these.d.hypers, these.q.hypers)
  
  plotPosteriorST(cat, cat.posts.plot, this.n=nsims, model = "spat-temp",
                  M0 = this.M0, 
                  hypers.list = this.hypers.list,
                  plot.path= paste0("plots/models/bayesianETAS/2020-experiments/st-model/", 
                                    this.plot.loc, "/posts-", cat.label, ".pdf"), 
                  priors.PB=this.PB.region, priors.JG=this.JG.region)
  # plot.path= paste0("plots/models/bayesianETAS/2020-experiments/st-model/", cat.string, "/multi-realns/", substr(cat.label, 1, 5), "/posts-", cat.label, ".pdf"))
  plotTraceST(cat.posts.plot, 
              plot.path = paste0("plots/models/bayesianETAS/2020-experiments/st-model/", 
                                 this.plot.loc, "/trace-", cat.label, ".pdf"))
  plotScatterMatrix(cat.posts.plot, plot.title = plot.title,
                    plot.path = paste0("plots/models/bayesianETAS/2020-experiments/st-model/", 
                                       this.plot.loc, "/smat-", cat.label, ".pdf"))
  # plot.path= paste0("plots/models/bayesianETAS/2020-experiments/st-model/", cat.string, "/multi-realns/", substr(cat.label, 1, 5), "/smat-", cat.label, ".pdf"))
  
}

# runSTBETASPlotsPNWPriors <- function(cat, cat.label, nsims, these.init.vals, this.M0,
#                                      this.spat.area = intl.t.area, this.isAlphaFixed=T,
#                                      cat.Gevals.mat, this.maxT, this.plot.loc,
#                                      plot.title, these.mu.hyperparams, these.K.hyperparams, these.alpha.hyperparams, 
#                                      these.c.hyperparams, these.p.hyperparams, these.d.hyperparams, these.q.hyperparams,
#                                      this.priors.PB) {
#   # cat=cat2D.nomod[cat2D.nomod$target,]; isCat2=T; isTarget=T; cat.label="cat2D-test"
#   # cat.Gevals.mat=cat2D.Gevals.mat; this.S.target=c(250,750)
#   # true.vals=cat2D.true; these.init.vals = abs(perturbInits(cat2D.true))
#   # plot.title = "Cat2D Test"
#   # nsims=5+burnin; maxT=15000
#   # 
#   # cat=intl.t; isCat2=F; isTarget=T; cat.label="pnw"
#   # cat.Gevals.mat=cat2D.Gevals; this.S.target=c(250,750)
#   # true.vals=cat2D.true; these.init.vals = abs(perturbInits(cat2D.true))
#   # plot.title = "Cat3B No mod, Target, Init Vals Near True"
#   # nsims=5+burnin; maxT=15000
#   # print("weeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee")
#   
#   cat.posts <- runSTBETASPrior(cat, init.vals=these.init.vals, 
#                                Geval.mat=cat.Gevals.mat, this.M0=M0,
#                                spatial.area=this.spat.area, sims=nsims, 
#                                maxT=this.maxT, isAlphaFixed = this.isAlphaFixed, 
#                                mu.hyperparams=these.mu.hyperparams,
#                                K.hyperparams=these.K.hyperparams,alpha.hyperparams=these.alpha.hyperparams,
#                                c.hyperparams=these.c.hyperparams,p.hyperparams=these.p.hyperparams,
#                                d.hyperparams=these.d.hyperparams,q.hyperparams=these.q.hyperparams)
#   write.csv(cat.posts, 
#             paste0("data/bayesianETAS-output/2020-experiments/st-model/pnw/posts-", cat.label, ".csv"), row.names=F)
#   cat.posts.plot <- cat.posts[-c(1:burnin),]
#   cat.posts.plot$mu <- cat.posts.plot$mu * this.spat.area
#   plotPosteriorST(cat, cat.posts.plot, this.n=nsims, model = "spat-temp",
#                   this.M0 = 2.0, priors.PB = this.priors.PB,
#                   plot.path= paste0("plots/models/bayesianETAS/2020-experiments/st-model/", 
#                                     this.plot.loc, "/posts-", cat.label, ".pdf"))
#   # plot.path= paste0("plots/models/bayesianETAS/2020-experiments/st-model/", cat.string, "/multi-realns/", substr(cat.label, 1, 5), "/posts-", cat.label, ".pdf"))
#   plotTraceST(cat.posts.plot, 
#               plot.path = paste0("plots/models/bayesianETAS/2020-experiments/st-model/", 
#                                  this.plot.loc, "/trace-", cat.label, ".pdf"))
#   plotScatterMatrix(cat.posts.plot, plot.title = plot.title,
#                     plot.path = paste0("plots/models/bayesianETAS/2020-experiments/st-model/", 
#                                        this.plot.loc, "/smat-", cat.label, ".pdf"))
#   # plot.path= paste0("plots/models/bayesianETAS/2020-experiments/st-model/", cat.string, "/multi-realns/", substr(cat.label, 1, 5), "/smat-", cat.label, ".pdf"))
#   
# }
# 
# runTempBETAS <- function(cat, init.vals){
#   # cat <- this.cat3A
#   # init.vals <- cat3A.true
#   mu.cat <- init.vals[1]; K.cat <- init.vals[2]; alpha.cat <- init.vals[3]; 
#   cp.cat <- init.vals[4]; p.cat <- init.vals[5];
#   out.cat.target.temp <- estimateETASBranchingInteractionMSList(ts=as.double(cat$time[cat$target]),
#                                                                 marks=as.double(cat$mag[cat$target]), 
#                                                                 # branching=as.integer(cat$reindex.branching[cat$target]),
#                                                                 branching=rep(0, length(cat$time[cat$target])),
#                                                                 maxT = as.double(maxT), 
#                                                                 M0=as.double(2.5), sims=as.integer(sims), numMCMCSamples=as.integer(numMCMCSamples),
#                                                                 mu=as.double(mu.cat)/(500^2), logK=as.double(log(K.cat)),
#                                                                 alpha=as.double(alpha.cat), c=as.double(cp.cat), p=as.double(p.cat),
#                                                                 mus=mus, logKs=logKs, alphas=alphas,
#                                                                 cs=cs, ps=ps)
#   Sys.time()
#   posts.cat.target.temp <- as.data.frame(cbind(out.cat.target.temp$nbkgd, out.cat.target.temp$mu, 
#                                                exp(out.cat.target.temp$logK), out.cat.target.temp$alpha,
#                                                out.cat.target.temp$c, out.cat.target.temp$p))
#   posts.cat.target.temp  <- posts.cat.target.temp[-c(1:burnin, nrow(posts.cat.target.temp)),]
#   
#   names(posts.cat.target.temp) <- c("nbkgd","mu", "K", "alpha", "c", "p")
#   # posts.cat.target.temp$mu <- posts.cat.target.temp$mu*(500^2)
#   # posts.cat.target.temp$mu <- posts.cat.target.temp$mu*500/((500^2)/(1000^2))
#   
#   return(posts.cat.target.temp)
#   
# }
