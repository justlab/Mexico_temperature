This is a repository for a model to predict temperature in an area around Mexico City, and the code to assemble a dataset for it.

We use several sources of data:

- Files for these networks are automatically downloaded by ``stations.R``.

  - SIMAT: The Sistema de Monitoreo Atmosférico de la Ciudad de México (SIMAT), which includes the Red de Meteorología y Radiación Solar (REDMET) and the Red de Depósito Atmosférico (REDDA). Website: http://www.aire.cdmx.gob.mx
  - UNAM: The Programa de Estaciones Meteorológicas del Bachillerato Universitario (PEMBU) at the Universidad Nacional Autónoma de México (UNAM). Website: https://www.ruoa.unam.mx/pembu

- Files from one family of networks have to be obtained in person, or by mailing a USB stick or hard drive. But information about ESIMEs and EMAs stations (specifically, their longitudes and latitudes) is automatically downloaded.

  - SMN: The Servicio Meteorológico Nacional México (SMN), which includes an apparently unnamed network of observatories (contact person: Adolfo Portocarrero Reséndiz; adolfo.portocarrero@conagua.gob.mx; (55) 2636-4600; Av. Observatorio 192, Col. Observatorio, Del. Miguel Hidalgo. C.P. 11860, México). It also includes the Estación Sinóptica Meteorológica (ESIMEs) network and the Estación Meteorológica Automática (EMAs) network (contact person: Lic. Moisés Espinosa Cárdenas; moises.espinosa@conagua.gob.mx; 01-(55)-26-36-46-00 extension 3484).

      - Since R's ``readxl`` package has difficulty with the EMAs Excel files, I use LibreOffice to mass-convert them to CSV with the following command: ``time ( ls | xargs --delimiter '\n' --max-args 250 --max-procs 10 soffice --headless --convert-to csv --outdir /data-belle/Mexico_temperature/stations/smn-emas-csv )``. This took about 1 hour 45 min on Belle.

- Files from Weather Underground can be downloaded or obtained from archives of previous downloads.

License
============================================================

This program is copyright 2018, 2019 Kodi B. Arfer and Iván Gutiérrez-Avila.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the `GNU General Public License`_ for more details.

.. _`GNU General Public License`: http://www.gnu.org/licenses/
