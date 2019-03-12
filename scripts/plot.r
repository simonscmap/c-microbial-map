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
  
  out_dir = normalizePath(out_dir)
  
  if (!dir.exists(out_dir)) {
    printf("Creating outdir '%s'\n", out_dir)
    dir.create(out_dir)
  }

  d <- read.csv(file_name)
  
  tmp <- split(d, d$latitude)
  names(tmp) <- paste0('stn', seq_along(tmp))
  
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
      oceSetData(ctd, 'eASV_Relative_Abundance', x$Relative_Abundance[o])
  })
  
  sec <- as.section(stns) #Create a section from the CTD data
  sg <-
    sectionGrid(sec, p = seq(20, 200, 5)) #For better interpolation (though be cautious with these plots!)
  
  base = tools::file_path_sans_ext(basename(file_name))
  outname <-
    file.path(out_dir, paste0(base, '-eASV-plot-%02d.png'))
  
  invisible(png(
    outname,
    width = img_width,
    height = img_height,
    res = 300,
    units = 'in'
  ))
  
  abun <-
    unlist(lapply(stns, function(x)
      x[['eASV_Relative_Abundance']]))
  
  #plot(sec, showstations=TRUE, ztype="contour", showBottom=FALSE)
  
  ###Dot plot with interpolated temperature up top, only colour proportional to the abundance
  par(mfrow = c(2, 1))
  plot(
    sg,
    which = 'temperature',
    xtype = "latitude",
    showstations = TRUE,
    showBottom = FALSE,
    ztype = 'image',
    zcol = oceColorsTemperature,
    cex = 2,
    pch = 16
  )
  cm <- colormap(abun, col = oceColorsViridis)
  #drawPalette(colormap=cm)
  plot(
    sec,
    which = 'eASV_Relative_Abundance',
    xtype = "latitude",
    showstations = TRUE,
    grid = TRUE,
    showBottom = FALSE,
    ztype = 'points',
    zcol = oceColorsViridis,
    cex = 2,
    pch = 16
  )
  
  #Change up the sizes
  #cex <- max(log10(abun))/log10(abun) #Normalized sizes
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
  
  ###Dot plot with interpolated temperature up top, color and the size proportional to the abundance
  par(mfrow = c(2, 1))
  
  plot(
    sg,
    which = 'temperature',
    xtype = "latitude",
    showstations = TRUE,
    showBottom = FALSE,
    ztype = 'image',
    zcol = oceColorsTemperature,
    cex = 2,
    pch = 16
  )
  
  plot(
    sec,
    which = 'eASV_Relative_Abundance',
    showBottom = FALSE,
    ztype = 'points',
    zcol = oceColorsViridis,
    cex = 0
  )
  
  # setup the plot
  points(sec[['distance']],
         sec[['pressure']],
         pch = 19,
         cex = 3 * cex,
         col = cm$zcol)
  
  ###All interpolated (sketchy!)
  par(mfrow = c(2, 1))
  
  plot(
    sg,
    which = 'temperature',
    xtype = "latitude",
    showstations = TRUE,
    showBottom = FALSE,
    ztype = 'image',
    zcol = oceColorsTemperature,
    cex = 2,
    pch = 16
  )
  
  plot(
    sg,
    which = 'eASV_Relative_Abundance',
    xtype = "latitude",
    showstations = TRUE,
    grid = TRUE,
    showBottom = FALSE,
    ztype = 'image',
    zcol = oceColorsViridis,
    cex = 2,
    pch = 16
  )
  
  invisible(dev.off())
  
  print('Done.')
}

main()