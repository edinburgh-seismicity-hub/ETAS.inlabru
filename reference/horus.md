# HORUS Ita Catalogue

The HOmogenized instRUmental Seismic catalog (HORUS) of Italy from 1960
to present, limited to limited to 1960-2019.

## Usage

``` r
horus
```

## Format

Original file format:

- `Year`::

  Origin Time (OT) year

- `Mo`::

  OT month

- `Da`::

  OT day

- `Ho`::

  OT hour

- `Mi`:

  OT minute

- `Se`:

  OT seconds and fractions

- `Lat`:

  epicenter N latitude (in degrees)

- `Lon`:

  epicenter E longitude (in degrees)

- `Depth`:

  hypocenter depth (in km)

- `Mw`:

  true or proxy moment magnitude

- `sigMw`:

  moment magnitude uncertainty

- `Geo.Ita`:

  `"*"` indicates that the epicenter is within the Italian mainland
  territory, otherwise `" "`

- `Geo.CPTI15`:

  `"*"` indicates that the epicenter is within the spatial window of the
  CPTI15 catalog (Rovida et al., 2020, Bull Earth Eng, doi:
  10.1007/s10518-020-00818-y)

- `Ev..type`:

  `"x"` indicates that the event is not an earthquake (e.g. explosion,
  eruption, landslide, ...) (only since May 1st 2012)

- `Iside.n.`:

  ISIDe id number (only since April 16th 2005)

ETAS.inlabru format:

- `time_string`::

  Data-time of the event, in ISO 8601 format;
  `format = "%Y-%m-%dT%H:%M:%OS"` for use with
  [`as.POSIXct()`](https://rdrr.io/r/base/as.POSIXlt.html)

- `lon`::

  Original `Lon`

- `lat`::

  Original `Lat`

- `depth`:

  Original `Depth`

- `M`::

  Original `Mw`

## Source

<http://horus.bo.ingv.it/>

## Details

The original entire HORUS catalog is provided as a single tab separated
ascii file: in `HORUS_Ita_Catalog.txt` a random second decimal digit is
added to ML and Md of the INGV bulletin (ISIDe) before computing Mw
proxies (Lolli et al., Seism. Res. Lett., 91, 3208-3222, doi:
10.1785/0220200148), in `HORUS_Ita_Catalog_o.txt` the original ML and Md
are used to compute Mw proxies. `ETAS.inlabru` includes a reformatted
version of the data from `HORUS_Ita_Catalog.txt` as a `data.frame`,
limited to 1960-2019.

The data are provided as they are, no express or implied warranty is
given.

## References

Barbara Lolli(1), Daniele Randazzo(1), Gianfranco Vannucci(1) and Paolo
Gasperini (2),(1) (2020). The Homogenized Instrumental Seismic Catalog
(HORUS) of Italy from 1960 to Present, Seismol. Res. Lett, doi:
10.1785/0220200148.

\(1\) Istituto Nazionale di Geofisica e Vulcanologia, Sezione di Bologna
(2) Dipartimento di Fisica e Astronomia, Università di Bologna
