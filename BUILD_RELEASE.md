# Notes on Building Releases

## Running Unit Tests

In order to run the unit tests manually do the following,

1. In R set the working directory to the `package` folder
in your local Git repository.
2. Make sure you have the R `devtools` package installed.
2. Execute the function `devtools::test()`.

## Building a CRAN release:

In order to build a release follow these steps:

 1. In your local repository go to the folder just above the `package` folder and do,
	```
	$ R CMD build package
	$ R CMD check --as-cran popsom_xyz.tar.gz
	```
	where `xyz` is the version number given in the DESCRIPTION file.  This builds a folder called `popsom.Rcheck`.

2. If the check passes, submit tarball `popsom_xyz.tar.gz` to `https://cran.r-project.org/submit.html`

**NOTE**: the `popsom` folder in the `popsom.Rcheck` folder is
a valid package. It can be loaded into R as follows:

	R> install.packages("<path to popsom.Rcheck/popsom>",repos=NULL)

## Building/Loading the Shareable Fortran Library

1. Go to the  `src` folder of the cloned git popsom repo and type the following command
```
$ R CMD SHLIB vsom.f90
```
This generates a .so and .o file.
1. Set the working directory of the R terminal to the location of vsom.so
1. Load `vsom.so` as a dynamic object using `dyn.load()`
1.  Verify that the dynamic loadable object is loaded using `is.loaded()`
1.  Load `map-utils.R` using `source()`. You will probably have to do
    some manual adjustments to the R code to load it outside of a package
    context.

Briefly:
```
R CMD SHLIB vsom.f90
dyn.unload("vsom.so")
dyn.load("vsom.so")
is.loaded("vsom.so")
```
