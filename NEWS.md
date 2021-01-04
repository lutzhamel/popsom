### POPSOM Release 5.0

- completely rearchitected package
- new, streamlined S3 based API
- popsom now supports two models:
  1. a self-organizing map model,
  2. a centroid based clustering model
- quality measures available for both models
- cleaned up map visualization
- extremely fast training algorithm based on ideas from tensor algebra

### Release 5.1

- Something got rattled with the S3 interface in 4.x.  It no longer works the way it did in release 3.x.  I took the S3 interface out because I want the package to work with both 3.x and 4.x installations.  
- Implemented a 'summary' function for map objects.
