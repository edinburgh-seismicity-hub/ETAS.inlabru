#' @name horus
#' @title HORUS Ita Catalogue
#' @docType data
#' @description The HOmogenized instRUmental Seismic catalog (HORUS) of Italy from 1960 to present
#'
#' @usage data(horus, package = "ETAS.inlabru")
#' @source <http://horus.bo.ingv.it/>
#'
#' @references
#' Barbara Lolli(1), Daniele Randazzo(1), Gianfranco Vannucci(1) and Paolo Gasperini (2),(1)
#' (2020). The Homogenized Instrumental Seismic Catalog (HORUS) of Italy from 1960 to Present,
#' Seismol. Res. Lett, doi: 10.1785/0220200148.
#'
#' (1) Istituto Nazionale di Geofisica e Vulcanologia, Sezione di Bologna
#' (2) Dipartimento di Fisica e Astronomia, Università di Bologna
#'
#' @details
#' The original entire HORUS catalog is provided as a single tab separated ascii file:
#'   in `HORUS_Ita_Catalog.txt` a random second decimal digit is added to ML and Md of the INGV
#' bulletin (ISIDe) before computing Mw proxies (Lolli et al., Seism. Res. Lett., 91,
#' 3208-3222, doi: 10.1785/0220200148), in `HORUS_Ita_Catalog_o.txt` the original ML and Md
#' are used to compute Mw proxies. `ETAS.inlabru` includes a reformatted version
#' of the data from `HORUS_Ita_Catalog.txt` as a `data.frame`.
#'
#' The data are provided as they are, no express or implied warranty is given.
#'
#' @format Original file format:
#' \describe{
#'   \item{`Year`:}{Origin Time (OT) year}
#'   \item{`Mo`:}{OT month}
#'   \item{`Da`:}{OT day}
#'   \item{`Ho`:}{OT hour}
#'   \item{`Mi`}{OT minute}
#'   \item{`Se`}{OT seconds and fractions}
#'   \item{`Lat`}{epicenter N latitude (in degrees)}
#'   \item{`Lon`}{epicenter E longitude (in degrees)}
#'   \item{`Depth`}{hypocenter depth (in km)}
#'   \item{`Mw`}{true or proxy moment magnitude}
#'   \item{`sigMw`}{moment magnitude uncertainty}
#'   \item{`Geo.Ita`}{`"*"` indicates that the epicenter is within the Italian
#'     mainland territory, otherwise `" "`}
#'   \item{`Geo.CPTI15`}{`"*"` indicates that the epicenter is within the spatial
#'     window of the CPTI15 catalog (Rovida et al., 2020, Bull Earth Eng, doi: 10.1007/s10518-020-00818-y)}
#'   \item{`Ev..type`}{`"x"` indicates that the event is not an earthquake
#'     (e.g. explosion, eruption, landslide, ...) (only since May 1st 2012)}
#'   \item{`Iside.n.`}{ISIDe id number (only since April 16th 2005)}
#' }
#' ETAS.inlabru format:
#' \describe{
#'   \item{`lon`:}{Original `Lon`}
#'   \item{`lat`:}{Original `Lat`}
#'   \item{`M`:}{Original `Mw`}
#'   \item{`time_String`:}{Combination of original `Year`, `Mo`, `Da`, `Mi`, and `Se`}
#'   \item{`depth`}{Original `Depth`}
#'   \item{`catalog_id`}{???}
#'   \item{`event_id`}{???}
#' }
NULL
