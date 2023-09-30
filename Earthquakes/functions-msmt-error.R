sampleCircle <- function(n, rad, center){
  # center <- current.eqk[, c("x", "y")]
  # rad <- next.gen.dists
  # n=current.eqk$num.trig
  
  center.x <- center$x ; center.y <- center$y
  u <- 2*pi*runif(n) # Random angle
  cbind(x=(rad*cos(u))+center.x, y=rad*sin(u)+center.y)
}

simCatsMsmtError <- function(this.cat, this.seed=NA, output.path = NA, output.path.ppm = NA, jitter.mags=T) {
  # this.seed <- 1616
  # this.cat <- intl.aux
  # jitter.mags <- T

  # Remove the x y that previously was in the cat (in case these had been recentered)
  this.cat$x <- this.cat$y <- NULL
  
  this.cat.xy <-  latlong2grid(this.cat[, c("lon", "lat")])
  this.cat <- cbind(this.cat, this.cat.xy)
  
  this.cat.x.range <- range(this.cat$x)
  this.cat.y.range <- range(this.cat$y)
  # Convert herr from km to degree y
  # this.cat$herr.deg <- 1/(this.cat$herr*111)
  
  if(!is.na(this.seed)) set.seed(this.seed)
  
  # Create new catalog with simulated ys, xs and mags
  # Horizontal error is actually the largest of the major axes of the error ellipsoid reported by the seismic network
  cat.sim <- this.cat
  # new.ys <- rnorm(n=nrow(this.cat), mean=this.cat$y, sd=this.cat$herr)
  # new.xs <- rnorm(n=nrow(this.cat), mean=this.cat$x, sd=this.cat$herr)
  new.locs <- as.data.frame(sampleCircle(n=nrow(this.cat), rad=this.cat$herr, center=data.frame("x"=this.cat$x, "y"=this.cat$y)))
  new.locs[is.na(this.cat$herr),] <- this.cat[is.na(this.cat$herr), c("x", "y")]
  # new.xs[is.na(new.locs)] <- this.cat$x[is.na(new.locs)]
  
  
  new.mags <- rnorm(n=nrow(this.cat), mean=this.cat$mag, sd=this.cat$merr)
  new.mags[is.na(new.mags)] <- this.cat$mag[is.na(new.mags)]
  cat.sim$x <- new.locs[,1]
  cat.sim$y <- new.locs[,2]

  if(jitter.mags) {
    cat.sim$mag <- round(new.mags, 1)
  }
  
  # Make an identifier if the y/x/mag is outside my spatio-mag experimental window
  cat.sim$out.of.win <- cat.sim$x > this.cat.x.range[2] | cat.sim$x < this.cat.x.range[1] | 
    cat.sim$y > this.cat.y.range[2] | cat.sim$y < this.cat.y.range[1] | 
    cat.sim$mag < min(this.cat$mag)
  
  # Sanity checks
  stopifnot(any(!is.na(cat.sim[, c("x", "y", "mag")])))
  
  # Make new lat/lon from these x/y
  cat.sim$lon <- cat.sim$lat <- NULL
  cat.sim.lat.lon <- grid2latlong(cat.sim[, c("x", "y")]) 
  names(cat.sim.lat.lon) <- c("lon", "lat")
  cat.sim <- cbind(cat.sim, cat.sim.lat.lon)
  
  # Recenter the x/y
  cat.sim$x <- (cat.sim$x - min(cat.sim$x)) + 1
  cat.sim$y <- (cat.sim$y - min(cat.sim$y)) + 1
  
  
  if(!is.na(output.path)) {
    write.csv(cat.sim, output.path, row.names=F)
    
    # Get rid of all the eqks outside of spatio-magnitude window before making HIST PPM catalog
    # cat.sim.in.window <- cat.sim[!cat.sim$out.of.win,]
    # stopifnot(nrow(cat.sim.in.window) == nrow(cat.sim)-sum(cat.sim$out.of.win))
    # cat.sim.ppm <- makeHISTPPMcatalog(cat.sim.in.window, start.year="1970")
    # write.table(cat.sim.ppm, output.path.ppm, sep="\t", 
    #             row.names=F, col.names = F)
  }
  
  return(cat.sim)
}

simCatSTBETASPNW <- function(init.cat, region.label, crustal=F, swarms=NA, target=F, fix.seed, random.seed, file.label){
  # init.cat <- intl.aux
  # region.label <- "north"
  # swarms <- "noCswarms"
  # target <- T
  # random.seed <- 2
  # file.label <- "Nt-noCswarms-Nt-cat2initsSS-fixalpha"
  
  j.init.cat <- simCatsMsmtError(this.cat=init.cat, jitter.mags=T, this.seed=fix.seed+random.seed)
  if(target){
    if(region.label == "north"){
      j.this.cat <- j.init.cat %>% filter(lat >= 45, lat <= 49, lon >= -125, lon <= -116.5, mag >= 2, 
                                          is.na(remove), date >= as.Date("1985-01-01"),
                                          potential_duplicates_20210412 == F)
    } else if(region.label == "south"){
      j.this.cat <- j.init.cat %>% filter(lat >= 42, lat < 45, lon >= -125, lon <= -116.5, mag >= 2, 
                                          is.na(remove), date >= as.Date("2004-01-01"),
                                          potential_duplicates_20210412 == F)
    } else if(region.label == "south-comp1"){
      j.this.cat <- j.init.cat %>% filter(lat >= 42, lat < 45, lon >= -125, lon <= -116.5, mag >= 2.3, 
                                          is.na(remove), date >= as.Date("2004-01-01"),
                                          potential_duplicates_20210412 == F)
    } else if(region.label == "south-comp2"){
      j.this.cat <- j.init.cat %>% filter(lat >= 42, lat < 45, lon >= -125, lon <= -116.5, mag >= 2.5, 
                                          is.na(remove), date >= as.Date("2004-01-01"),
                                          potential_duplicates_20210412 == F)
    }
  }
  
  file.name <- paste0("data/processed-data/pnw-cats-msmt-error/", region.label, "/jLM", random.seed, "-", file.label, ".csv")
  if(swarms == "noCswarms"){
    j.this.cat <- j.this.cat %>% filter(!scheme2020 %in% c(1))
  } else if(swarms == "noMCswarms"){
    j.this.cat <- j.this.cat %>% filter(!scheme2020 %in% c(1, 2))
  }
  
  if(crustal){
    j.this.cat <- j.this.cat %>% filter(Depth2ModS > 10| is.na(Depth2ModS))
  }
  
  j.this.cat$x <- j.this.cat$x - min(j.this.cat$x)
  j.this.cat$y <- j.this.cat$y - min(j.this.cat$y)
  j.this.cat <- j.this.cat[order(j.this.cat$time),]
  j.this.cat$time <- j.this.cat$time - min(j.this.cat$time)
  write.csv(j.this.cat, file.name)
  
  return(j.this.cat)
} 
