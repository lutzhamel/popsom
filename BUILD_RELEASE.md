## Building a CRAN release:

You can find the popsom github repo here:

	https://github.com/lutzhamel/popsom

If you don't have a copy of the popsom repo clone it:

	$ git clone https://github.com/lutzhamel/popsom.git

 otherwise do a pull to make sure you have the latest version:

 	$ git pull

In order to build a release do the following:

 1. Go to the folder just above the local popsom repo and do:
	```
	$ R CMD build popsom
	$ R CMD check --as-cran popsom_xyz.tar.gz
	```
	where `xyz` is the version number given in the DESCRIPTION file.  This builds a folder called `popsom.Rcheck`.

2. Submit tarball `popsom_xyz.tar.gz` to `https://cran.r-project.org/submit.html`

**NOTE**: the `popsom` folder in the `popsom.Rcheck` folder is
a valid package. It can be loaded into R as follows:

	R> install.packages("<path to popsom.Rcheck/popsom>",repos=NULL)

**NOTE**: to build a shareable Fortran library for R:
1. Go to the  `src` folder of the cloned git popsom repo.
1. Input `R CMD SHLIB vsom.f90` into the Linux terminal, which outputs a .so and .o file.
1. Set the working directory of the R terminal to the location of vsom.so
1. Load `vsom.so` as a dynamic object using dyn.load()
1.  Verify that the dynamic loadable object is loaded using `is.loaded()`
1.  Load `map-utils.R` using `source()`

Briefly:

	R CMD SHLIB vsom.f90

	dyn.unload("vsom.so")
	dyn.load("vsom.so")
