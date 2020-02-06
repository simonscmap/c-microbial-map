#!/usr/bin/env Rscript

# Authors: Jesse McNichol <mcnichol@alum.mit.edu>
#          Ken Youens-Clark <kyclark@email.arizona.edu>

suppressMessages(library("oce"))
suppressMessages(library("optparse"))
suppressMessages(library("R.utils"))

# -----------------------------------------------------------
getArgs = function() {
  args = commandArgs(trailingOnly = TRUE)

  option_list = list(
    make_option(
      c("-f", "--file"),
      default = "",
      type = "character",
      help = "Input file",
      metavar = "character"
    ),
    make_option(
      c("-o", "--outdir"),
      default = file.path(getwd(), 'plots'),
      type = "character",
      help = "Output directory",
      metavar = "character"
    ),
    make_option(
      c("-t", "--title"),
      default = "",
      type = "character",
      help = "Title",
      metavar = "character"
    ),
    make_option(
      c("-l", "--legend"),
      default = "",
      type = "character",
      help = "Legend",
      metavar = "character"
    ),
    make_option(
      c("-w", "--width"),
      default = 7,
      type = "integer",
      help = "Image width",
      metavar = "int"
    ),
    make_option(
      c("-H", "--height"),
      default = 7,
      type = "integer",
      help = "Image height",
      metavar = "int"
    )
  )
  opt_parser = OptionParser(option_list = option_list)
  return(parse_args(opt_parser))
}

# -----------------------------------------------------------
main = function() {
  opts = getArgs()
  file_name  = opts$file
  out_dir    = opts$outdir
  img_width  = opts$width
  img_height = opts$height
  plot_title = opts$title
  plot_legend = opts$legend

  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }

  out_dir = normalizePath(out_dir)

  d <- read.csv(file_name)

  aStations = seq_along(split(d, d$longitude)) #create array of length of unique longitude  values
  keeporder <- factor(d$longitude, levels=unique(d$longitude))
  tmp <- split(d, keeporder)
  names(tmp) <- paste0('stn', aStations)

  #Import the stations as a CTD object
  stns <- lapply(tmp, function(x) {
    o <- order(x$depth) # the depths are not always in order
    ctd <-
      as.ctd(
        x$salinity[o],
        x$temperature[o],
        swPressure(x$depth[o]),
        longitude = x$longitude,
        latitude = x$latitude
      )
    ctd <- oceSetData(ctd, 'depth', x$depth[o])
    ctd <-
      oceSetData(ctd, 'eASV_Relative_Abundance', x$relative_abundance[o])
  })

  sec <- lon360(as.section(stns)) #Create a section from the CTD data, transforming to 360 degree notation since GP13 crosses date line
  maxDepth <- max(sec[['depth']])
  rescaledDepth <- maxDepth * 1.25

  sg <-
    sectionGrid(sec, p = seq(20, 200, 5)) #For better interpolation (though be cautious with these plots!)

  base = tools::file_path_sans_ext(basename(file_name))
  outname <-
    file.path(out_dir, paste0(base, '-eASV-plot-%02d.png'))

  invisible(png(
    outname,
    width = img_width,
    height = img_height,
    res = 600,
    units = 'in'
  ))

  abun <-
    unlist(lapply(stns, function(x)
      x[['eASV_Relative_Abundance']]))

  ###Dot plot with interpolated temperature up top, only colour proportional to the abundance
  par(mfrow = c(3, 1), oma=c(2,2,2,2))

  plot(
    sg,
    which = 'temperature',
    xtype = "longitude",
    showstations = TRUE,
    showBottom = FALSE,
    ztype = 'image',
    zcol = oceColorsTemperature,
    cex = 2,
    pch = 16
  )

  if (length(plot_title) > 0) {
    mtext(plot_title, 3, line=0, outer=TRUE, cex=1.2)
  }

  cm <- colormap(abun, col = oceColorsViridis)
  #drawPalette(colormap=cm)
  plot(
    sec,
    which = 'eASV_Relative_Abundance',
    xtype = "longitude",
    showstations = TRUE,
    grid = TRUE,
    showBottom = FALSE,
    ztype = 'points',
    zcol = oceColorsViridis,
    cex = 2,
    pch = 16,
    ylim = c(rescaledDepth,0)
  )

  if (length(plot_legend) > 0) {
    mtext(plot_legend, 1, line=1, outer=TRUE, cex=0.8)
  }

  plot(
    sec,
    which = 99,
    showstations=TRUE,
    showStart=TRUE
  )


  #Change up the sizes
  cex <- abun / max(abun) #Normalized sizes

  ###T-S diagram
  par(mfrow = c(1, 1))
  cm <- colormap(abun, col = oceColorsViridis)
  drawPalette(colormap = cm)
  plotTS(
    sec,
    pch = 19,
    col = cm$zcol,
    mar = c(3.5, 3.5, 2, 4),
    cex = 3 * cex
  )

  if (length(plot_title) > 0) {
    mtext(plot_title, 3, line=0, outer=TRUE, cex=1.2)
  }

  if (length(plot_legend) > 0) {
    mtext(plot_legend, 1, line=1, outer=TRUE, cex=0.8)
  }

  ###Dot plot with interpolated temperature up top, color and the size proportional to the abundance
  par(mfrow = c(3, 1))

  plot(
    sg,
    which = 'temperature',
    xtype = "distance",
    showstations = TRUE,
    showBottom = FALSE,
    ztype = 'image',
    zcol = oceColorsTemperature,
    cex = 2,
    pch = 16
  )

  if (length(plot_title) > 0) {
    mtext(plot_title, 3, line=0, outer=TRUE, cex=1.2)
  }

  plot(
    sec,
    which = 'eASV_Relative_Abundance',
    showBottom = FALSE,
    ztype = 'points',
    zcol = oceColorsViridis,
    cex = 0,
    ylim = c(rescaledDepth,0)
  )

  # setup the plot
  points(sec[['distance']],
         sec[['pressure']],
         pch = 19,
         cex = 3 * cex,
         col = cm$zcol)

#  data(coastlineWorld)
#  cl <- coastlineCut(coastlineWorld, -180)
#  mapPlot(cl, col="lightgray",
#        projection="+proj=lcc +lat_0=-30 +lon_0=165 +lat_1=-20 +lat_2=-40",
#        longitudelim=c(150,200), latitudelim=c(-45,-10))
#  mapPoints(sec, cex=0.5)


   plot(
     sec,
     which = 99,
     showstations=TRUE,
     showStart=TRUE
   )

   if (length(plot_legend) > 0) {
     mtext(plot_legend, 1, line=1, outer=TRUE, cex=0.8)
   }

  ###All interpolated (sketchy!)
  par(mfrow = c(3, 1))

  plot(
    sg,
    which = 'temperature',
    xtype = "longitude",
    showstations = TRUE,
    showBottom = FALSE,
    ztype = 'image',
    zcol = oceColorsTemperature,
    cex = 2,
    pch = 16
  )

  if (length(plot_title) > 0) {
    mtext(plot_title, 3, line=0, outer=TRUE, cex=1.2)
  }

  plot(
    sg,
    which = 'eASV_Relative_Abundance',
    xtype = "longitude",
    showstations = TRUE,
    grid = TRUE,
    showBottom = FALSE,
    ztype = 'image',
    zcol = oceColorsViridis,
    cex = 2,
    pch = 16
  )

  plot(
    sec,
    which = 99,
    showstations=TRUE,
    showStart=TRUE
  )

  if (length(plot_legend) > 0) {
    mtext(plot_legend, 1, line=1, outer=TRUE, cex=0.8)
  }

  invisible(dev.off())

  printf("Done, see outdir '%s'\n", out_dir)
}

main()
