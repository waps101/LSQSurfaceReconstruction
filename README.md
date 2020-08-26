# [Least squares surface reconstruction on arbitrary domains](https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123670528.pdf)

 [Dizhong Zhu](https://www.cs.york.ac.uk/cvpr/member/dizhong/) and [William A. P. Smith](https://www-users.cs.york.ac.uk/wsmith)
 
 #### [ECCV2020]

<br/>

## Abstract

Almost universally in computer vision, when surface derivatives are required, they are computed using only first order accurate finite difference approximations. We propose a new method for computing numerical derivatives based on 2D Savitzky-Golay filters and K-nearest neighbour kernels. The resulting derivative matrices can be used for least squares surface reconstruction over arbitrary (even disconnected) domains in the presence of large noise and allowing for higher order polynomial local surface approximations. They are useful for a range of tasks including normal-from-depth (i.e. surface differentiation), height-from-normals (i.e. surface integration) and shape-from-x. We show how to write both orthographic or perspective height-from-normals as a linear least squares problem using the same formulation and avoiding a nonlinear change of variables in the perspective case. We demonstrate improved performance relative to state-of-the-art across these tasks on both synthetic and real data and make available an open source implementation of our method.

## Matlab implementation

The key function is gradMatricesSavGol which takes as input a 2D binary mask and returns derivative and smoothing matrices to compute derivatives for each foreground pixel using a Savitzky-Golay filter of specified size and order. The function is used as follows:
```matlab
    [Dx,Dy,Serr,S] = gradMatricesSavGol(mask,nOrder,nSize,verbose)
```
where mask is the rows x cols binary foreground mask, nOrder is a scalar specifying the order of polynomial surface to use (typical values are 2..4), nSize is the width of the filter when used as a square domain (should be odd, typical values 3..7), verbose is a flag to specify whether some information is displayed. Dx and Dy are sparse matrices of size npix x npix where npix is the number of true values in mask. They compute horiztonal and vertical derivatives when multiplied by the vectorised version of the foreground function values, e.g. Dx * z will compute the horizontal derivatives of the function values in z, stored as a vector. To put the values back into the matrix you could do something like:
```matlab
    p = zeros(size(mask));
    p(mask) = Dx * z;
```
The other two outputs are related to smoothing and are of the same size as the derivative matrices. S computes the smoothed function. Serr computes the difference between the smoothed and original function.

As an example application of the derivative matrices, we also include a function for doing surface integration (aka height/depth-from-gradient or height/depth-from-normals). This is called as follows:
```matlab
    [z] = surfaceFromGradient(Dx,Dy,p,q,mask,options)
```
Dx and Dy are the derivative matrices computed above and mask is the binary foreground mask. p and q are the input noisy derivatives. options allows a number of options to be passed in. A typical setup for perspective integration might be:
```matlab
    options.structure = 'stacked';
    options.precond = 'none';
    options.solver = 'direct';
    options.regularisation='smoothness';
    options.lambda = 2;
    options.cameramodel = 'perspective';
    options.cx = cx;
    options.cy = cy;
    options.fx = fx;
    options.fy = fy;
```
See the function for the other options that can be used.

## Citation

If you use the model or the code in your research, please cite the following paper:

Dizhong Zhu and William A. P. Smith. "Least squares surface reconstruction on arbitrary domains". In Proc. of the European Conference on Computer Vision (ECCV), 2020.
[https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123670528.pdf](https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123670528.pdf)

Bibtex:

    @inproceedings{smith2020morphable,
      title={Least squares surface reconstruction on arbitrary domains},
      author={Zhu, Dizhong and Smith, William A. P.},
      booktitle={Proc. of the European Conference on Computer Vision (ECCV)},
      year={2020}
    }
