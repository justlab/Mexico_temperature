suppressPackageStartupMessages(
   {library(data.table)
    library(sf)
    library(ggplot2)})

source("common.R")

l = local(
   {source("modeling.R")
    list(master.grid, ground, stations, predict.temps, per.mrow.population, all.agebs.year)})
master.grid = l[[1]]
ground = l[[2]]
stations = l[[3]]
predict.temps = l[[4]]
per.mrow.population = l[[5]]
all.agebs.year = l[[6]]

base.size = 14

temp.quantiles.map = function(the.year)
   {d = local(
      {p = predict.temps(the.year, "pred.area")
       p = rbind(
           p[, .(mrow, temp = pred.ground.temp.lo, kind = 1)],
           p[, .(mrow, temp = pred.ground.temp.hi, kind = 2)])
       cbind(p,
           master.grid[p$mrow, .(x_sinu, y_sinu)])})
    d = d[, by = .(x_sinu, y_sinu, kind), .(temp = quantile(temp, .95))]
    message("Range: ", floor(min(d$temp)), " to ", ceiling(max(d$temp)))
    ggplot(d) +
        geom_raster(aes(x_sinu, y_sinu, fill = temp)) +
        facet_grid(. ~ kind) +
        scale_fill_distiller(
            name = "Temperature (°C)",
            palette = "Spectral",
            limits = c(0, 40),
            breaks = seq(0, 40, by = 5),
            guide = guide_colorbar(nbin = 500)) +
        geom_vline(data = data.frame(kind = 1, x = max(d$x_sinu) + 1),
            aes(xintercept = x)) +
        theme_void(base_size = base.size) +
        theme(
           strip.background = element_blank(),
           strip.text.x = element_blank(),
           panel.spacing = unit(-7, "mm"),
           legend.position = "bottom", legend.key.width = unit(20, "mm")) +
        coord_equal()}

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
        theme_void() +
        coord_equal()}

area.map = function(years = NULL)
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
            data = master.grid[unique(stations[
                (if (is.null(years)) T else (stn %in%
                    ground[year(date) %in% years, unique(stn)])),
                mrow])]) +
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

mexico.context.map = function()
    ggmap::ggmap(get_stamenmap(zoom = 5,
            c(left = -123, right = -86, bottom = 10, top = 34),
            maptype = "terrain")) +
        with(study.area(), annotate("rect",
            xmin = left, xmax = right,
            ymin = bottom, ymax = top,
            color = "red", fill = NA)) +
        theme_void()

pop.map = function(pop.col, thresholds.tempC = NULL)
   {d = merge(by = "mrow", master.grid, per.mrow.population(pop.col))
    if (!is.null(thresholds.tempC))
       {d = merge(d, by = "mrow", all.x = T,
            predict.temps(all.agebs.year, "pred.area")[, keyby = mrow, .(
                extreme.days = c(
                    sum(pred.ground.temp.lo <= thresholds.tempC[1]),
                    sum(pred.ground.temp.hi >= thresholds.tempC[2])),
                kind = c(
                    sprintf("≤ %d °C", thresholds.tempC[1]),
                    sprintf("≥ %d °C", thresholds.tempC[2])))])
        stopifnot(!anyNA(d$extreme.days))
        d[, val := extreme.days * pop]
        print(d[, keyby = kind, .("total person-days" = sum(val))])}
    else
        setnames(d, "pop", "val")

    p = ggplot() +
        geom_raster(aes(x_sinu, y_sinu, fill = val),
            data = d) +
        theme_void(base_size = base.size) +
        coord_equal()
    if (is.null(thresholds.tempC)) p = p +
        scale_fill_distiller(
            name = "Population",
            palette = "RdPu", direction = 1,
            guide = guide_colorbar(nbin = 500),
            breaks = c(0, 500, 1000, 1500, 2000, 2700),
            limits = c(0, 2700),
            labels = scales::comma)
    else p = p +
        scale_fill_distiller(
            name = "Exposure\n(person-days)",
            palette = "RdPu", direction = 1,
            guide = guide_colorbar(nbin = 500),
            breaks = seq(0, 140e3, 20e3),
            limits = c(0, 140e3),
            labels = function(x) ifelse(x == 0, "0",
                paste0(round(x/1000), "k"))) +
        facet_grid(. ~ kind) +
        geom_vline(data = data.frame(kind = d$kind[1], x = max(d$x_sinu) + 1),
            aes(xintercept = x)) +
        theme(
           strip.background = element_blank(),
           panel.spacing = unit(-7, "mm"),
           legend.position = "bottom", legend.key.width = unit(20, "mm"))

    p}

time.series.plot = function()
   {stns = c("Mexico City" = 8, Morelos = 24)
    years = c(2010L, 2018L)
    mon = 6

    # Get cross-validated predictions.
    d = rbindlist(lapply(years, function(y)
        (run.cv(y, "ground.temp.mean")
            [stn %in% unname(stns),
                c(.SD, .(date = lubridate::make_date(y, 1, 1) - 1 + yday))]
            [month(date) == mon,
                .(stn, year, mday = mday(date), pred, ground.temp)])))
    # Verify that we have a prediction for every day at `stns` in
    # `years`.
    d[, by = .(stn, year), stopifnot(.N ==
        lubridate::days_in_month(lubridate::make_date(year, mon, 1)))]
    # Reshape the data.
    d = melt(d, id.vars = c("stn", "year", "mday"),
        variable.name = "type", value.name = "temp")
    d[, type := factor(ifelse(type == "pred", "Predicted", "Observed"))]
    d[, stn := factor(names(stns)[match(stn, stns)])]

    ggplot(d, aes(mday, temp, color = type, group = type)) +
        geom_point() +
        geom_line() +
        facet_grid(stn ~ year) +
        scale_x_continuous(name = "Day of June",
            breaks = c(1, 10, 20, 30)) +
        scale_y_continuous(name = "Temperature (°C)",
            breaks = seq(15, 30, by = 5),
            limits = c(15, 30),
            expand = c(0, 0)) +
        labs(color = "") +
        theme_bw() +
        theme(
           axis.text = element_text(color = "black"),
           panel.grid.major.x = element_blank(),
           panel.grid.minor = element_blank())}

pred.error.plot = function()
   {the.year = 2018L

    d = run.cv(the.year, "ground.temp.mean")[!is.na(pred)]
    d[, season := months2seasons[month(
        lubridate::make_date(the.year, 1, 1) - 1 + yday)]]

    ggplot(d) +
        geom_density(aes(pred - ground.temp), bw = .1) +
        facet_grid(season ~ .) +
        xlab("Prediction error") +
        scale_y_continuous(limits = c(0, .4)) +
        coord_cartesian(xlim = c(-10, 5), expand = F) +
        theme_bw() +
        theme(
           axis.text = element_text(color = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.title.y = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank())}
