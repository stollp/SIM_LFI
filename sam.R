# Function to  summarize measurements for a single terrestrial sample plots at xs/ys

#  xs/ys: x/y-coordinate of sample plot center
#      r: radius of sample plot = 12.62 m -> 500 m2

sam <- function(xs = 1000, ys = 1000, r = 12.62, Ps = P, OUTPUT = TRUE) {

  xl <- xs - r; xu <- xs + r; yl <- ys - r; yu <- ys + r
  
  Ps <- Ps[Ps$X >= xl & Ps$X <= xu & Ps$Y >= yl & Ps$Y <= yu, ]
  
dist <- sqrt((Ps$X - xs)^2 + (Ps$Y - ys)^2)  # distance to sample plot center

  Ps <- data.frame(Ps, dist); Ps <- Ps[dist <= r, ] 

  if (dim(Ps)[1] == 0 & OUTPUT) cat("\nno trees in sample plot\n\n")
  
#  if (dim(Ps)[1]  > 0 & OUTPUT) {print(Ps)}
  
#  cat(dim(Ps)[1] * 20, mean(Ps[, 'B'], sum(Ps[, 'B']) * 20), "\n")
  
  return(list(dim(Ps)[1] * 20, sum(Ps[, 'B']) * 20, sum(Ps[, 'B'] * 20))) # scale to hectare
  
}