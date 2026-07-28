pphotspot
=====================

https://cran.r-project.org/package=pphotspot

The `R` package `pphotspot` provides implements the detection of hotspots
on a linear network as proposed by Mrkvička et al. (2025). 
The package utilizes the `R` package `GET`(see Myllymäki and Mrkvička, 2024).

## The development version

The github repository holds a copy of the current development version of the contributed R package `pphotspot`.

This development version is as or more recent than the official release of `GET` on the Comprehensive R Archive Network (CRAN) at https://cran.r-project.org/package=GET

## Where is the official release?

For the most recent official release of `pphotspot`, see https://cran.r-project.org/package=pphotspot

## Installation

### Installing the official release

To install the official release of `pphotspot` from CRAN, start `R` and type

```R
install.packages('pphotspot')
```

### Installing the development version

The easiest way to install the `pphotspot` library from github is through the `remotes` package. Start `R` and type:

```R
require(remotes)
install_github('myllym/pphotspot')
```
If you do not have the R library `remotes` installed, install it first by running

```R
install.packages("remotes")
```

After installation, in order to start using `pphotspot`, load it to R and
see the main help page, which describes the functions of the library:
```R
require(pphotspot)
help('pphotspot-package')
```

If you want to have also vignettes working, you should also install packages from the 'suggests' field,
have MiKTeX or tinytex on your computer, and install the library with
```R
install_github('myllym/pphotspot', build_vignettes = TRUE)
```

## Vignettes

The package contains a vignette. The vignette shows the package usage.
It is available by starting `R` and typing
```R
library("pphotspot")
vignette("HotSpots")
```

The vignette is also available at the package webpage https://cran.r-project.org/package=pphotspot

## Branches

Currently the is only the main branch called `master`.

## References

To cite GET in publications use

Mrkvička, T., Kraft, S., Blažek, V., Myllymäki, M., Konopa, M. (2025)
Where are the Real Hotspots of Road Crashes? An Innovative Approach to the Road Crash Hotspot Detection by Filtering Out the Effect of Covariates. 
Available at SSRN: http://dx.doi.org/10.2139/ssrn.5337003

Myllymäki, M. and Mrkvička, T. (2024).
GET: Global envelopes in R. Journal of Statistical Software 111(3), 1-40. doi: 10.18637/jss.v111.i03 https://doi.org/10.18637/jss.v111.i03
