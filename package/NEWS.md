### Release 6.0

- Renamed the functions in the interface to avoid collisions with S3 functions within the R environment. We know that this is another renaming of the POPSOM interface and we apologize for any inconvenience.  We expect that the interface is now stable for the foreseeable future.
- New features:
    * The `map.minimal` object.  This is an object that only contains the trained neurons and nothing else. This is an appropriate model when POPSOM is used as a preprocessing step and no other model information is needed.  Note that `map.minimal` objects cannot be processed by any of the other functions in the POPSOM interface.
    * The `map.convergence` function provides details about the underlying convergence characteristics.
- Bugfixes
  * Most importantly the artificial limit of a minimum of 50 instances in the training
    data has been removed.

### Release 5.2

Reworked the description of the package in order to reflect the capabilities
of the package better.

### Release 5.1

- Something got rattled with the S3 interface in 4.x.  It no longer works the way it did in release 3.x.  I took the S3 interface out because I want the package to work with both 3.x and 4.x installations.  
- Implemented a 'summary' function for map objects.

### Release 5.0

- completely rearchitected package
- new, streamlined S3 based API
- popsom now supports two models:
  1. a self-organizing map model,
  2. a centroid based clustering model
- quality measures available for both models
- cleaned up map visualization
- extremely fast training algorithm based on ideas from tensor algebra
