myHeatmap <-
#fastest way to project a matrix to a binary bitmap
function(x,ord,xlab="",ylab="",main="My Heatmap",
                      col=heat.colors(5), ...){
    op <- par(mar=c(3,0,2,0)+0.1)
    on.exit(par(op))
    nc <- NCOL(x)
    nr <- NROW(x)
    labCol <- names(x)

    x <- t(x[ord,])
    image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
        c(0, nr), axes = FALSE, xlab=xlab, ylab=ylab, main=main,
        col=col,...)

    axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0)
    axis(2, 1L:nr, labels = NA, las = 2, line = -0.5, tick = 0)
}
