suppressPackageStartupMessages(
   {library(data.table)
    library(sf)
    library(ggplot2)})

l = local(
   {source("modeling.R")
    list(master.grid, stations, pred.area, predict.temps, per.mrow.population, crs.satellite, all.agebs.year)})
master.grid = l[[1]]
stations = l[[2]]
pred.area = l[[3]]
predict.temps = l[[4]]
per.mrow.population = l[[5]]
crs.satellite = l[[6]]
all.agebs.year = l[[7]]

temp.map = function(the.year, temp.kind, agg, fill.args)
   {d = local(
      {p = predict.temps(the.year, "pred.area")
       cbind(
           p[, .(temp = get(paste0("pred.ground.temp.", temp.kind)))],
           master.grid[p$mrow, .(x_sinu, y_sinu)])})
    d = d[, by = .(x_sinu, y_sinu), .(temp = agg(temp))]
    message("Range: ", floor(min(d$temp)), " to ", ceiling(max(d$temp)))
    ggplot(d) +
        geom_raster(aes(x_sinu, y_sinu, fill = temp)) +
        do.call(scale_fill_distiller, c(fill.args, list(
            palette = "Spectral",
            limits = range(fill.args$breaks),
            guide = guide_colorbar(nbin = 500)))) +
        theme_void()}

change.map = function(years1, years2)
   {d = local(
      {f = function(ys)
           rbindlist(lapply(ys, function(y) predict.temps(y, "pred.area")))[,
               keyby = mrow, .(temp = mean(pred.ground.temp.mean))]
       p1 = f(years1)
       p2 = f(years2)
       cbind(
           change = p2$temp - p1$temp,
           master.grid[p1$mrow, .(x_sinu, y_sinu)])})
    message("Range: ", round(min(d$change), 2), " to ", round(max(d$change), 2))
    ggplot(d) +
        geom_raster(aes(x_sinu, y_sinu, fill = change)) +
        scale_fill_gradient2(
            low = scales::muted("blue"),
            mid = "#ffffcc",
            high = scales::muted("red"),
            limits = c(-2.7, 1.6),
            breaks = c(-2.7, -2, -1, 0, 1, 1.6),
            guide = guide_colorbar(nbin = 500)) +
        theme_void()}

area.map = function()
   {xr = range(master.grid$x_sinu)
    yr = range(master.grid$y_sinu)
    xd = 5000
    yd = 5000
    ggplot() +
        geom_raster(aes(x_sinu, y_sinu), fill = "gray95",
            data = master.grid) +
        geom_sf(data = pred.area(),
            fill = "white", color = "black", size = .2) +
        geom_point(aes(x_sinu, y_sinu), color = "red",
            size = .1,
            data = master.grid[unique(stations$mrow)]) +
        coord_sf(crs = crs.satellite, expand = F,
            xlim = c(xr[1] - xd, xr[2] + xd),
            ylim = c(yr[1] - yd, yr[2] + yd)) +
        scale_x_continuous(name = "",
            breaks = -100 : -98) +
        scale_y_continuous(name = "",
            breaks = c(18.5, 19, 19.5, 20, 20.5)) +
        theme_bw() +
        theme(
            panel.grid.major = element_line(color = "transparent"))}

pop.map = function(pop.col, threshold.tempC = NULL, hotter = T, xlims = NULL, ylims = NULL)
   {d = merge(
        master.grid[if (is.null(xlims) | is.null(ylims)) T else
           (xlims[1] <= lon & lon <= xlims[2] &
            ylims[1] <= lat & lat <= ylims[2])],
        per.mrow.population(pop.col),
        by = "mrow")
    if (!is.null(threshold.tempC))
       {d = merge(d,
            predict.temps(all.agebs.year, "pred.area")[, keyby = mrow,
                .(hot.days = sum(if (hotter)
                    pred.ground.temp.hi >= threshold.tempC else
                    pred.ground.temp.lo <= threshold.tempC))],
            by = "mrow",
            all.x = T)
        stopifnot(!anyNA(d$hot.days))
        d[, val := hot.days * pop]}
    else
        setnames(d, "pop", "val")

    ggplot() +
        geom_raster(aes(x_sinu, y_sinu, fill = val),
            data = d) +
        scale_fill_distiller(
            palette = "Spectral",
            guide = guide_colorbar(nbin = 500)) +
        theme_void() +
        coord_equal()}
