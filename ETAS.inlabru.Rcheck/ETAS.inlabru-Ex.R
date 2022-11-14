pkgname <- "ETAS.inlabru"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "ETAS.inlabru-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('ETAS.inlabru')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("create.input.list.temporal.noCatalogue")
### * create.input.list.temporal.noCatalogue

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: create.input.list.temporal.noCatalogue
### Title: Function to create a default input file for the ETAS Hawkes
###   temporal model where no catalogue is specified in the input file
### Aliases: create.input.list.temporal.noCatalogue

### ** Examples

# HOW DO WE REFERENCE A FILE IN THE data DIRECTORY?
#create.input.list.temporal.noCatalogue('data/user_input_synthetic_noCatalog.txt')



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("create.input.list.temporal.noCatalogue", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("generate.temporal.ETAS.synthetic")
### * generate.temporal.ETAS.synthetic

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: generate.temporal.ETAS.synthetic
### Title: Generates a sythetic catalogue using the ETAS model
### Aliases: generate.temporal.ETAS.synthetic

### ** Examples

## EXAMPLE 1: Generate a 1000 day synthetic ETAS catalogue

generate.temporal.ETAS.synthetic( theta=data.frame(mu=0.1, K=0.089, alpha=2.29, c=0.11, p=1.08), beta.p=log(10), M0=2.5, T1=0, T2=1000 )


## EXAMPLE 2: To generate a 1000 day catalogue including a M6.7 event on day 500

Ht <- data.frame(ts=c(500), magnitudes=c(6.7))
generate.temporal.ETAS.synthetic( theta=data.frame(mu=0.1, K=0.089, alpha=2.29, c=0.11, p=1.08), beta.p=log(10), M0=2.5, T1=0, T2=1000, Ht=Ht )



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("generate.temporal.ETAS.synthetic", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sample.GR.magnitudes")
### * sample.GR.magnitudes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sample.GR.magnitudes
### Title: Return a sample of magnitudes drawn from the GR distribution
### Aliases: sample.GR.magnitudes

### ** Examples

sample.GR.magnitudes(n=100, beta.p=log(10), M0=2.5)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sample.GR.magnitudes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sample.temporal.ETAS.generation")
### * sample.temporal.ETAS.generation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sample.temporal.ETAS.generation
### Title: Take all previous parent events from 'Ht=data.frame[ts,
###   magnitudes]' and generates their daughters events using the ETAS
###   model
### Aliases: sample.temporal.ETAS.generation

### ** Examples

# The parents are specified in Ht
Ht <- data.frame(ts=c(500), magnitudes=c(6.7))
sample.temporal.ETAS.generation( theta=data.frame(mu=0.1, K=0.089, alpha=2.29, c=0.11, p=1.08), beta.p=log(10), M0=2.5, T1=0, T2=1000, Ht=Ht )



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sample.temporal.ETAS.generation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("time.grid")
### * time.grid

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: time.grid
### Title: Generate a set of time bins for a specific event and return.
### Aliases: time.grid

### ** Examples

## EXAMPLE 1
events <- data.frame( ts=c(0,1 , 3 ), idx.p=c(1,2,3) )
T2 <- 20
N.exp <- 8
delta.t <- 0.1
coef.t <- 1
time.grid(events, coef.t, delta.t, T2, displaygrid = FALSE, N.exp)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("time.grid", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
