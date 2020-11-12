In this project, "A spatiotemporal reconstruction of daily ambient temperature using satellite data in the Megalopolis of Central Mexico from 2003–2018", we built a model to predict the mean, maximum, and minimum temperature on each day at each square in a 1-km grid for an area around Mexico City.

Raw and processed data, including predictions, as well as a research notebook, can be found on Zenodo at http://doi.org/10.5281/zenodo.3362523

The bulk of the data-processing and data-analysis code can be found on GitHub at https://github.com/justlab/Mexico_temperature . Other code is embedded inline in the manuscript and notebook.

Instructions
============================================================

Getting the temperature predictions
------------------------------------------------------------

To examine or use our temperature predictions without running any of our code, take a look at the `HDF5 <http://portal.hdfgroup.org/display/HDF5/Introduction+to+HDF5>`_ files ``predictions_*.h5``. There's one file per year. Each file has a three-dimensional array ``data`` with an attribute ``dimensions`` naming the dimensions (location, time, and variable) and a group ``dimension_labels`` naming each index of each dimension. The temperatures are in degrees Celsius, and the dates indicate 24-hour spans of UTC-06:00. ``mrow`` values are row indices of the master grid, which can be found in ``master_grid.h5``. The original coordinates of the grid (in which it is, in fact, a regular grid) are ``x_sinu`` and ``y_sinu``, which are in the coordinate reference system ``crs.satellite``, defined in ``common.R`` as ``"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"``, which is the `MODIS global sinusoidal projection <https://modis-land.gsfc.nasa.gov/MODLAND_grid.html>`_.

Here's how you could plot the mean temperatures for 5 July 2012 in R

.. code-block:: R

    library(hdf5r)
    library(ggplot2)

    h5 = H5File$new("predictions_2012.h5", mode = "r")
    preds = h5[["data"]][,,]
    dimnames(preds) = sapply(simplify = F,
        h5attr(h5[["data"]], "dimensions"),
        function(k) h5[["dimension_labels"]][[k]][])
    preds = preds[, "2012-07-05", "pred.ground.temp.mean"]
    h5$close_all()

    h5 = H5File$new("master_grid.h5", mode = "r")
    g = h5[["master_grid"]][]
    h5$close_all()
    g$pred = preds[as.character(g$mrow)]

    ggplot(g[!is.na(g$pred),]) +
        geom_raster(aes(x_sinu, y_sinu, fill = pred))

Getting the ground-station observations
------------------------------------------------------------

We put a lot of effort into cleaning and unifying daily observations of air temperature (and wind speed, precipitation, and air pressure) at ground weather stations from several different sources. We provide all the original raw data on Zenodo, but if you'd like to use the processed observations without running any of the processing code, use the `JSON <https://www.json.org>`_ file ``ground.json.gz``, which has information about the observations as well as the stations they come from. To open this file in R, run ``jsonlite::fromJSON(gzfile("ground.json.gz"))``. Again, dates are in UTC-06:00.

Running the code
------------------------------------------------------------

To reproduce our results, you have the option of starting from scratch, that is, generating ``ground.json.gz`` from the raw station data, or of using the cleaned observations in ``ground.json.gz`` and modeling temperature from there. Either way, you'll need to:

1. Install any libraries required by ``library(...)`` calls in the file you're using.
2. Set the environment variable ``JUSTLAB_MEXICO_TEMPERATURE_DATA_ROOT`` to the directory containing the data from Zenodo.

To generate ``ground.json.gz``, source ``stations.R``. You'll need to download and uncompress the ``geography`` and ``stations`` data from Zenodo. You can then call ``save.ground()``.

To perform cross-validation and generate predictions, source ``modeling.R``. You'll need ``geography`` (uncompressed) and ``ground.json.gz`` (left compressed) from Zenodo. You can then do cross-validation and summarize the results with a call like ``summarize.cv.results(run.cv(2012L, "ground.temp.mean"))`` or get all the predictions for a year with a call like ``predict.temps(2012L, "pred.area")``.

Notes
============================================================

Since R's ``readxl`` package has difficulty with the EMAs Excel files, I used LibreOffice to mass-convert them to CSV with the following command: ``time ( ls | xargs --delimiter '\n' --max-args 250 soffice --headless --convert-to csv --outdir ../../smn-emas-csv/2018 )`` (and likewise ``smn-emas-csv/2019``). This can take a while.

License
============================================================

We assert no copyright over the data. The code is copyright 2018, 2019, 2020 Kodi B. Arfer, Iván Gutiérrez-Avila, and Johnathan Rush.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the `GNU General Public License`_ for more details.

.. _`GNU General Public License`: http://www.gnu.org/licenses/
