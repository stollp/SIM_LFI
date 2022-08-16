# Script to simulate population and sampling -> estimates

source("sam.R")

set.seed(13)

require( plotrix, quietly = TRUE)
require(EnvStats, quietly = TRUE)  

graphics.off()
# output                unknown Pop full census/ha
	   GRAPH <- TRUE; POP <- TRUE; FULL <- TRUE

		 n.s <- 16			# sample size, 49 yields 500 x 500  m grid
							#			   16          1 x   1 km
 	sampling <- 'Grid'
#  	sampling <-  'SRS'

		   N <- 1 * 10^6	# unknown number of trees

		   w <- h <- 10		# dimension (w * h [km]) of unknown population

B <- rlnorm(     n = N, mean = 2,       sd = 1)		# ba size distribution m²
B <- rlnormTrunc(n = N, meanlog = 0.01, sdlog = 0.4, min = .01, max = 10^6)/10

r <- 100 * sqrt(B/pi)								# radius [cm]
D <-   2 * r          								#    dbh [cm]
U <-  (D * pi)/100                          			#    gbh [ m]

X <- runif(n = N,  min = 0, max = w^4)/2   			# tree positions [m]
Y <- runif(n = N,  min = 0, max = h^4)/2   			# aggregated in llc

P <- data.frame( X, Y, r, D, U, B)         			#; print(head(P))
P[, 1:5] <- round(P[, 1:5], 1)

if (POP) { "on/off"

	cat("\nUnknown, infinite population\n\n")

	cat("          Total area [Km^2]:", sprintf("%8.0f",  w * h),        "      Forested area [Km^2]:", 
	                                    sprintf("%5.0f", (w * h)/4),                    "\n", sep = '')
	cat("            Number of trees:", sprintf("%8.0f",      N),                     "\n\n", sep = '')
	cat("         Number of trees/ha:", sprintf("%8.2f",      N/ 2500),                 "\n", sep = '')
	cat("    Mean basal area [m2/ha]:", sprintf("%8.2f",  sum(B)/2500),                 "\n", sep = '')
	cat("      Total basal area [m2]:", sprintf("%8.2f",  sum(B)/10^4), " (* 10^4 m2)", "\n", sep = '')

}

# count nr trees/ha & calc mean(ba/ha) and their variances

x <-    cut(P$X, breaks = 50)
y <-    cut(P$Y, breaks = 50)	# yields 2500 1 ha (100 * 100 m) tiles

T <- tapply(P$B, INDEX = list(x, y), FUN = length); 
M <- tapply(P$B, INDEX = list(x, y), FUN =    sum)	# schon auf Ha on funktion sam hochgerechnet

O <- (M * 2500)/10^4				# Total auf Forested erea hochgerechnet

if (length(T[is.na(T)]) > 0) cat("\n some NaN's replaced with 0\n")

T[is.na(T)] <- 0					# because some might be empty if N < 10^6
M[is.na(M)] <- 0
O[is.na(O)] <- 0

# Varianz: von jeder einzelnen Ha auf's Total hochrechnen (*2500), und daraus
#		   mean bzw. Varianz schaetzen.

if (FULL) { "on/off"

	cat("\nEstimates based on 2500 fully censused 1 ha tiles\n\n")

n.se   <-       sd(T)/sqrt(2500)		# abs. se's
b.se   <-       sd(M)/sqrt(2500)
t.se   <-       sd(O)/sqrt(2500)		

n.se.p <- 100/mean(T) * n.se 		# rel. se's
b.se.p <- 100/mean(M) * b.se
t.se.p <- 100/mean(O) * t.se		

# check whether CI covers unknown population
n.c <- mean(T) - 2 * n.se <     N /2500 &
       mean(T) + 2 * n.se >     N /2500
b.c <- mean(M) - 2 * b.se < sum(B)/2500 & 
       mean(M) + 2 * b.se > sum(B)/2500
t.c <- mean(O) - 2 * t.se < sum(B)	
       mean(O) + 2 * t.se > sum(B)

	cat("                                            95%-CI covers unknown population\n\n")

	cat("    Mean number of trees/ha:", sprintf("%8.2f", mean(T)),                  
	            "                 se:", sprintf("%6.2f",    n.se), 
	                              " (", sprintf("%4.1f", 	n.se.p), "%) ", sprintf("%5s", n.c), "\n", sep = '')
	cat("    Mean basel area [m2/ha]:", sprintf("%8.2f", mean(M)),                  
	            "                 se:", sprintf("%6.2f",    b.se),
	                              " (", sprintf("%4.1f",    b.se.p), "%) ", sprintf("%5s", b.c), "\n", sep = '')
	cat("      Total basal area [m2]:", sprintf("%8.2f", mean(O)), " (* 10^4 m2)", 
	                        "     se:", sprintf("%6.2f",    t.se), 
	                              " (", sprintf("%4.1f",    t.se.p), "%) ", sprintf("%5s", t.c), "\n", sep = '')

}

if (sampling == 'SRS') {	# simple random sampling

 X.S <- runif(n = n.s, min = 1000, max = 4000)
 Y.S <- runif(n = n.s, min = 1000, max = 4000)

}

if (sampling == 'Grid') {
# systematic sampling, 16 * 16
# X.S <- seq(from = 1000, to = 4000, by = 1000)
# Y.S <- seq(from = 1000, to = 4000, by = 1000); n.s <- length(X.S) * length(Y.S)

  X.S <- seq(from = 1000, to = 4000, by = (4000 - 1000)/(sqrt(n.s) - 1))
  Y.S <- seq(from = 1000, to = 4000, by = (4000 - 1000)/(sqrt(n.s) - 1)); n.s <- length(X.S) * length(Y.S)

# X.S <- X.S + 250   					# shift 'random' starting point of grid
# Y.S <- Y.S + 250   					# cf Magnussen et al. 2020

  X.S <- rep(X.S,        length(X.S))	# make grid
  Y.S <- rep(X.S, each = length(Y.S))

    g <- X.S[2] - X.S[1]				# distance between grid cells

}

if (sampling == 'SRS')  cat("\nEstimates based on", sprintf("%5.0f", n.s) , " terrestrial sample plots (simple random sampling)\n\n", sep = '')
if (sampling == 'Grid') cat("\nEstimates based on", sprintf("%5.0f", n.s) , " terrestrial sample plots (", sprintf("%3.2f", g/1000), " * ", sprintf("%3.2f", g/1000), " Km grid)\n\n", sep = '')

S <- data.frame(nr.trees = NaN, total.ba = NaN, Total.BA = NaN)			# containing the samples

for (i in 1: length(X.S)) {
  
    S[i, ] <- sam(xs = X.S[i], ys = Y.S[i], r = 12.62, Ps  = P, OUTPUT = FALSE)
    
  }

S[, 'Total.BA'] <- (S[, 'Total.BA'] * 2500)/10^4

n.se   <- sd(S[, 'nr.trees'])/sqrt(n.s)			# abs. se's
b.se   <- sd(S[, 'total.ba'])/sqrt(n.s) 
t.se   <- sd(S[, 'Total.BA'])/sqrt(n.s)

n.se.p <- 100/mean(S[, 'nr.trees']) * n.se 		# rel. se's
b.se.p <- 100/mean(S[, 'total.ba']) * b.se 		
t.se.p <- 100/mean(S[, 'Total.BA']) * t.se 		

# check whether CI covers unknown population
n.c <- mean(S[, 'nr.trees']) - 2 * n.se <     N /2500 & 
	   mean(S[, 'nr.trees']) + 2 * n.se >     N /2500
b.c <- mean(S[, 'total.ba']) - 2 * b.se < sum(B)/2500 & 
       mean(S[, 'total.ba']) + 2 * b.se > sum(B)/2500
t.c <- mean(S[, 'Total.BA']) - 2 * t.se < sum(B)/10^4 &	 
       mean(S[, 'Total.BA']) + 2 * t.se > sum(B)/10^4

cat("                                            95%-CI covers unknown population\n\n")
cat("   Mean number of trees/ha: ", sprintf("%8.2f", mean(S[, 'nr.trees'])),                        
           "                 se: ", sprintf("%5.2f", n.se), " (", sprintf("%4.1f", n.se.p), "%) ", sprintf("%5s", n.c), "\n", sep = '')
cat("   Mean basel area [m2/ha]: ", sprintf("%8.2f", mean(S[, 'total.ba'])),                        
           "                 se: ", sprintf("%5.2f", b.se), " (", sprintf("%4.1f", b.se.p), "%) ", sprintf("%5s", b.c), "\n", sep = '')
cat("     Total basal area [m2]: ", sprintf("%8.2f", mean(S[, 'Total.BA'])), " (* 10^4 m2)", 
                       "     se: ", sprintf("%5.2f", t.se), " (", sprintf("%4.1f", t.se.p), "%) ", sprintf("%5s", t.c), "\n", sep = '')

if (GRAPH) {
  
x11(width = 7, height = 7, pointsize = 16, title = 'test', xpos = 0, ypos = 0); par(bty = 'n', tck = 0.01)

plot(X, Y, pch = '.', xaxs = 'i', yaxs = 'i', main = expression(paste('100 ', Km^2)),
    xlab = 'X-Coordinate [Km]', xlim = c(0, 10000), ylim = c(0, 10000),
    ylab = 'Y-Coordinate [Km]', xaxt = 'n', yaxt = 'n', col = 'green')

  axis(1, at = c( 0, 100, 1000, 2000, 5000, 10000), lab = c(0, 0.1, 1, 2, 5, 10))
  axis(2, at = c( 0, 100, 1000, 2000, 5000, 10000), lab = c(0, 0.1, 1, 2, 5, 10))

  for (i in seq(from = 0, to = 10000, by = 1000)) abline(h = i, col = 'grey')
  for (i in seq(from = 0, to = 10000, by = 1000)) abline(v = i, col = 'grey')
 
  text(5000, 7000, 'non-forest', cex = 1, srt=-30)

  for (i in 1: length(X.S)) {
    
      points(X.S[i], Y.S[i], col = 'blue', cex = .5) # draw sample plot centers

    }

  }

# x11(width = 7, height = 7, pointsize = 16, title = '1 ha', xpos =   0, ypos = 700); par(bty = 'n', tck = 0.01); hist(B); 

if (FALSE) cat("\ndbg range [cm]", sprintf("%5.1f", range(D)), "   gbh range [m]", sprintf("%5.1f", range(U)), "\n")

if (GRAPH) {
# quartz()
# 700
x11(width = 7, height = 7, pointsize = 16, title = '1 ha', xpos =   0, ypos = 700); par(bty = 'n', tck = 0.01)

# Zoom
z <- 100; P <- P[P$X < z & P$Y < z, ]

plot(P$X, P$Y, type = 'n', xaxs = 'i', yaxs = 'i', main = expression(paste('1 Ha')),
     xlab = 'X-Coordinate [Km]', xlim = c(0, z), ylim = c(0, z),
     ylab = 'Y-Coordinate [Km]', xaxt = 'n', yaxt = 'n')

legend(-1, 55, col = c('red', 'blue', 'blue'), bty = 'n', cex = .5, lty = c(1, 1, 2),
       leg = c(expression(paste('Interpretation (50 * 50 m)')),
               expression(paste('Sample (500 ', m^2, ', r = 12.62 m)')),
               expression(paste('Sample (200 ', m^2, ', r =   7.98 m)'))))

#                                                            cm -> m
for (i in 1:dim(P)[1]) draw.circle(P$X[i], P$Y[i], radius = P$r[i]/100, col = 'green')

axis(1, at = c( 0, 50, 100, 1000, 2000, 5000, 10000), lab = c(0, '50 m', 0.1, 1, 2, 5, 10))
axis(2, at = c( 0, 50, 100, 1000, 2000, 5000, 10000), lab = c(0, '50 m', 0.1, 1, 2, 5, 10))

for (i in seq(from = 0, to = 10000, by = 1000)) abline(h = i, col = 'grey', xaxs = 'i', yaxs = 'i')
for (i in seq(from = 0, to = 10000, by = 1000)) abline(v = i, col = 'grey', xaxs = 'i', yaxs = 'i')

lines(c( 0, 50), c( 0,  0), lwd = 1, col = 'red') # interpretation area
lines(c( 0, 50), c(50, 50), lwd = 1, col = 'red')
lines(c( 0,  0), c( 0, 50), lwd = 1, col = 'red')
lines(c(50, 50), c( 0, 50), lwd = 1, col = 'red')

points(25,     25, pch = '+', col = 'blue', lwd = 1)
draw.circle(x = 25, y = 25, radius =  7.98, border = 'blue', lty = 2) # 200 m2 for trees >= 12 cm dbh < 36 cm
draw.circle(x = 25, y = 25, radius = 12.62, border = 'blue', lty = 1) # 500 m2 for trees >= 36 cm dbh
}
