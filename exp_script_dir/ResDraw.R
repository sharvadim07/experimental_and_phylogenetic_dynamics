options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
# args[1] - oct res file
# args[2] - param. variation config. file
# args[3] - file with coordinates of mode for the set of varied parameters
# args[4] - output file name

library (ggplot2)

dens_plot <- function (dat, file, x1, y1, x2, y2, x3, y3, x4, y4, mean_x, mean_y, sigma_x, sigma_y) {
    pdf (file)
    minx = min (dat [,1], x1, x2, x3)
    maxx = max (dat [,1], x1, x2, x3)

    miny = min (dat [,2], y1, y2, y3)
    maxy = max (dat [,2], y1, y2, y3)

    num_sample <- 100000
    pr_x <- rnorm (num_sample, mean=mean_x, sd=sigma_x)
    pr_y <- rnorm (num_sample, mean=mean_y, sd=sigma_y)    

    data_pr <- data.frame (pr_x, pr_y)
    colnames (data_pr) <- c ("x", "y")

    plot0 <- ggplot(data_pr, aes(x=x, y=y)) + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE) +
                   geom_point (aes (x=x1,y=y1), color="gray", size=5, shape=14) + 
                   scale_x_continuous(expand = c(0, 0), limits = c (mean_x - 3*sigma_x, mean_x + 3*sigma_x)) +
                   scale_y_continuous(expand = c(0, 0), limits = c (mean_y - 3*sigma_y, mean_y + 3*sigma_y)) +
                   scale_fill_distiller(palette="Spectral") +
                   xlab(conf [ind [j], 1]) + ylab(conf [ind [k], 1]) +
                   ggtitle("A priori probability density") + theme(plot.title = element_text(hjust = 0.5))
    print (plot0)


    plot1 <- ggplot(dat, aes(x=x, y=y)) + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE) +
                   geom_point (aes (x=x1,y=y1), color="gray", size=5, shape=14) + 
                   geom_point (aes (x=x2,y=y2), color="white", size=5, shape=9) + 
                   geom_point (aes (x=x3,y=y3), color="white", size=5, shape=8) + 
                   geom_point (aes (x=x4,y=y4), color="white", size=5, shape=0) + 
                   scale_x_continuous(expand = c(0, 0), limits = c (minx, maxx)) +
                   scale_y_continuous(expand = c(0, 0), limits = c (miny, maxy)) +
                   scale_fill_distiller(palette="Spectral") +
                   xlab(conf [ind [j], 1]) + ylab(conf [ind [k], 1]) +
                   ggtitle("A posteriori probability density (wide bounding box)") + theme(plot.title = element_text(hjust = 0.5))
    print (plot1)

    plot2 <- ggplot(dat, aes(x=x, y=y)) +  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
                   geom_point (aes (x=x1,y=y1), color="gray", size=5, shape=14) + 
                   geom_point (aes (x=x2,y=y2), color="white", size=5, shape=9) + 
                   geom_point (aes (x=x3,y=y3), color="white", size=5, shape=8) + 
                   geom_point (aes (x=x4,y=y4), color="white", size=5, shape=0) + 
                   scale_x_continuous(expand = c(0, 0), limits = c (minx, maxx)) +
                   scale_y_continuous(expand = c(0, 0), limits = c (miny, maxy)) +
                   scale_fill_distiller(palette="Spectral") +
                   xlab(conf [ind [j], 1]) + ylab(conf [ind [k], 1]) +
                   ggtitle("A posteriori probability density (wide bounding box)") + theme(plot.title = element_text(hjust = 0.5))

    print (plot2)

    xsorted <- sort (dat [,1])
    xlen <- length (xsorted)
    minx <- xsorted [as.integer (0.05 * xlen)]
    maxx <- xsorted [as.integer (0.95 * xlen)]

    ysorted <- sort (dat [,2])
    ylen <- length (ysorted)
    miny <- ysorted [as.integer (0.05 * ylen)]
    maxy <- ysorted [as.integer (0.95 * ylen)]

    plot3 <- ggplot(dat, aes(x=x, y=y)) + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE) +
                   geom_point (aes (x=x1,y=y1), color="gray", size=5, shape=14) + 
                   geom_point (aes (x=x2,y=y2), color="white", size=5, shape=9) + 
                   geom_point (aes (x=x3,y=y3), color="white", size=5, shape=8) + 
                   geom_point (aes (x=x4,y=y4), color="white", size=5, shape=0) + 
                   scale_x_continuous(expand = c(0, 0), limits = c (minx, maxx)) +
                   scale_y_continuous(expand = c(0, 0), limits = c (miny, maxy)) +
                   scale_fill_distiller(palette="Spectral") +
                   xlab(conf [ind [j], 1]) + ylab(conf [ind [k], 1]) +
                   ggtitle("A posteriori probability density (smaller bounding box)") + theme(plot.title = element_text(hjust = 0.5))
    print (plot3)

    plot4 <- ggplot(dat, aes(x=x, y=y)) +  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
                   geom_point (aes (x=x1,y=y1), color="gray", size=5, shape=14) + 
                   geom_point (aes (x=x2,y=y2), color="white", size=5, shape=9) + 
                   geom_point (aes (x=x3,y=y3), color="white", size=5, shape=8) + 
                   geom_point (aes (x=x4,y=y4), color="white", size=5, shape=0) + 
                   scale_x_continuous(expand = c(0, 0), limits = c (minx, maxx)) +
                   scale_y_continuous(expand = c(0, 0), limits = c (miny, maxy)) +
                   scale_fill_distiller(palette="Spectral") +
                   xlab(conf [ind [j], 1]) + ylab(conf [ind [k], 1]) +
                   ggtitle("A posteriori probability density (smaller bounding box)") + theme(plot.title = element_text(hjust = 0.5))
    print (plot4)


    tmp <- dev.off ()

    return ("OK")
}


dat_all <- read.table (args [1], header=FALSE)
conf <- read.table (args [2], header=FALSE, stringsAsFactors=FALSE)

pnum <- dim (conf) [1]
if (pnum != 16) {
    stop ("Error: invalid number of parameters in the param.conf. file\n")
}
for (j in 1:pnum) {
    qname <- paste("a", j, sep="")
    if (conf [j,1] != qname) {
        stop ("Error: invalid parameter name or order in the param.conf. file")
    }    
}

ind <- which (conf [,2] != "fixed")
num_ind <- length (ind)

mode_r <- as.vector (read.table (args [3], header=FALSE))
if (length (mode_r) != num_ind) {
    stop ("Error: inconsistent number of components in the configuration file and mode file\n")
}


conf [1,1] <- "log10 (Birth rate)"
conf [2,1] <- "log10 (Cmax)"
conf [3,1] <- "log10 (WT inf.rate)"
conf [4,1] <- "log10 (DI inf.rate)"
conf [5,1] <- "log10 (uninf.cell death rate)"
conf [6,1] <- "log10 (inf.cell death rate)"
conf [7,1] <- "log10 (WT binding rate)"
conf [8,1] <- "log10 (DI binding rate)"
conf [9,1] <- "log10 (WV production by a WC)"
conf [10,1] <- "log10 (DV production by a WDC)"
conf [11,1] <- "log10 (WV production by a WDC)"
conf [12,1] <- "log10 (WT inf.rate for RCs)"
conf [13,1] <- "log10 (DI inf.rate for RCs)"
conf [14,1] <- "log10 (Num. cells rec. signal after WT infection)"
conf [15,1] <- "log10 (Num. cells rec. signal after DI infection)"
conf [16,1] <- "log10 (Fraction of RCs in descendents of RCs)"


# generate point estimates of parameters
priors_r <- log10 (conf [, 3])

mean_r <- vector (mode="numeric", length=num_ind)
median_r <- vector (mode="numeric", length=num_ind)

for (j in 1:num_ind) {
    mean_r [j] <- mean (dat_all [, j])
    median_r [j] <- median (dat_all [, j])
}
    
res_r <- priors_r
res_r [ind] <- mode_r []
res_r <- as.numeric (res_r, mode="numeric")
cat (10^res_r, sep="\t", file=args[4], append=FALSE)
cat ("\n", file=args[4], append=TRUE)

res_r [ind] <- median_r
res_r <- as.numeric (res_r, mode="numeric")
cat (10^res_r, sep="\t", file=args[4], append=TRUE)
cat ("\n", file=args[4], append=TRUE)

res_r [ind] <- mean_r
res_r <- as.numeric (res_r, mode="numeric")
cat (10^res_r, sep="\t", file=args[4], append=TRUE)
cat ("\n", file=args[4], append=TRUE)

# Prepare density plots

if (num_ind < 2) {
    # Nothing to do
    q ()
}

labels_for_filenames <- c ("a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10", "a11", "a12", "a13", "a14", "a15", "a16")


for (j in 1:(num_ind-1)) {
    for (k in (j+1):num_ind) {
        fname <- paste (labels_for_filenames [ind [j]], "_", labels_for_filenames [ind [k]], ".pdf", sep="")

        dat <- dat_all [, c (j, k)]
        colnames (dat) <- c ("x", "y")

        x1 <- log10 (conf [ind [j], 3])
        y1 <- log10 (conf [ind [k], 3])

        x2 <- mean (dat [,1])
        y2 <- mean (dat [,2])

        x3 <- median (dat [,1])
        y3 <- median (dat [,2])

        x4 <- mode_r [j]
        y4 <- mode_r [k]

        mean_x <- x1
        mean_y <- y1
        sigma_x <- conf [ind [j], 4] / 3
        sigma_y <- conf [ind [k], 4] / 3

        tmp <- dens_plot (dat=dat, file=fname, x1=x1, y1=y1, x2=x2, y2=y2, x3=x3, y3=y3, x4=x4, y4=y4, mean_x=mean_x, mean_y=mean_y, sigma_x=sigma_x, sigma_y=sigma_y)
    }
}

