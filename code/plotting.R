suppressPackageStartupMessages(
   {library(data.table)
    library(sf)
    library(ggplot2)})

l = local(
   {source("modeling.R")
    list(master.grid, stations, pred.area, predict.temps)})
master.grid = l[[1]]
stations = l[[2]]
pred.area = l[[3]]
predict.temps = l[[4]]

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
   {ggplot() +
        geom_raster(aes(x_sinu, y_sinu), fill = "#aaaaaa",
            data = master.grid) +
        geom_sf(data = pred.area(),
            fill = "white", color = "black", size = .2) +
        geom_point(aes(x_sinu, y_sinu), color = "red",
            size = .1,
            data = master.grid[unique(stations$mrow)]) +
        coord_sf(crs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs",
            datum = NA) +
        theme_void()}