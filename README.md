## Fast 2D phase unwrapping implementation in MATLAB

Fast unwrapping 2D phase image using the algorithm given in:              

> M. A. Herraez, D. R. Burton, M. J. Lalor, and M. A. Gdeisat,          
  "Fast two-dimensional phase-unwrapping algorithm based on sorting by  
  reliability following a noncontinuous path", Applied Optics, Vol. 41,
  Issue 35, pp. 7437-7444 (2002).                                       

If using this code for publication, please kindly cite the following:
* M. A. Herraez, D. R. Burton, M. J. Lalor, and M. A. Gdeisat, *"Fast
two-dimensional phase-unwrapping algorithm based on sorting by reliability
following a noncontinuous path"*, Applied Optics, Vol. 41, Issue 35, pp. 7437-7444 (2002).
* M. F. Kasim, *"Fast 2D phase unwrapping implementation in MATLAB"*,
https://github.com/mfkasim91/unwrap_phase/ (2017).

### Getting started

* To get started and see the demo, simply type `test_unwrap_phase` in MATLAB
  (make sure that the file is included in MATLAB search path)
* The function to use is `unwrap_phase`
* Examples of using the code is given in the `test_unwrap_phase.m` file

### If this algorithm fails ...

Please try the [`unwrap_phase`](https://scikit-image.org/docs/stable/api/skimage.restoration.html?highlight=unwrap_phase#skimage.restoration.unwrap_phase) by scikit-image in Python. It has a better implementation.
