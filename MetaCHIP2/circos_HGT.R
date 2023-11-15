
# check.packages function: install and load multiple R packages.
# Check to see if packages are installed. Install them if they are not, then load them into the R session.
check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = 1)
}


# install packages if not 
packages<-c("optparse", "circlize")
invisible(suppressMessages(check.packages(packages)))

options(warn=-1)
option_list = list(
  make_option(c("-m", "--matrix"),  type="character", metavar="character", help="input matrix"),
  make_option(c("-s", "--fontsize"),type="double",    default=12,          help="font size, default: 12"),
  make_option(c("-p", "--plot"),    type="character", metavar="character", help="output plot"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

font_size = opt$fontsize/10

mat = read.table(opt$matrix, header = TRUE)

pdf(opt$plot, width=20, height=20, pointsize=12)
grid.col = c(A = 'brown1', B = 'lawngreen', C = 'mediumorchid', D = 'mediumslateblue', E = 'royalblue', F = 'sandybrown')
par(mar = rep(0,4), cex = font_size)
label_order = sort(union(rownames(mat), colnames(mat)))
chordDiagram(t(mat), order = label_order, annotationTrack = "grid", preAllocateTracks = 1, grid.col = grid.col)

# rotate label
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = font_size)
  #circos.axis(h = "top", labels.cex = font_size, sector.index = sector.name, track.index = 2, major.tick = TRUE, labels = TRUE)  # labels.cex control the size of scale
  circos.axis(h = "top", labels.cex = par("cex"), sector.index = sector.name, track.index = 2, major.tick = TRUE, labels = TRUE)  # labels.cex control the size of scale
}, bg.border = NA)

#circos.track(circos.text, facing = "clockwise")
invisible(dev.off())
rm(list=ls())
