__author__ = "Uwe Schmitt uschmitt@mineway.de, wrapped by V.A. Sole - ESRF"
__doc__ = """
This module is a simple wrapper to the py_nnma module of Uwe Schmitt (uschmitt@mineway.de)
in order to integrate it into PyMca. What follows is the documentation of py_nnma

py_nnma:  python modules for nonnegative matrix approximation (NNMA)

(c) 2009 Uwe Schmitt, uschmitt@mineway.de

NNMA minimizes  dist(Y, A X)

       where:  Y >= 0,  m x n
               A >= 0,  m x k
               X >= 0,  n x k

               k < min(m,n)

     dist(A,B) can be || A - B ||_fro 
                   or   KL(A,B)


This moudule provides the following functions:

    NMF, NMFKL, SNMF, RRI, ALS, GDCLS, GDCLS_L1, FNMAI, FNMAI_SPARSE,
    NNSC and FastHALS

The common parameters when calling such a function are:

    input:

            Y           --   the matrix for decomposition, maybe dense 
                             from numpy or sparse from scipy.sparse 
                             package

            k           --   number of componnets to estimate

            Astart 
            Xstart      --   matrices to start iterations. Maybe None
                             for using random start matrices.

            eps         --   termination swell value

            maxcount    --   max number of iterations to be performed

            verbose     --   if False: produce no output durint interations
                             if integer: give all 'verbose' itetations some
                             output about current state of iterations

    output:

            A, X        --   result matrices of algorithm 
            
            obj         --   value of objective function of last iteration

            count       --   number of iterations done

            converged   --   flag: indicates if iterations stoped within 
                             max number of iterations

The following extra parameters exist depending on algorithm:

    RRI      :  damping parameter 'psi' (default: 1e-12)

    SNMF     :  sparsity parameter 'sparse_par' (default: 0)

    ALS      :  regularization parameter 'regul' for stabilizing iterations
                (default value 0). needed if objective value jitters.
 
    GCDLS    :  'regul' for l2-smoothness of X (default 0)

    GDCLS_L1 :  'regul' for l1-smoothness of X (default 0)

    FNMAI    :  'stabil' for stabilizing algorithm (default value 1e-12)
                'alpha'  for stepsize  (default value 0.1)
                'tau'    for number of inner iterations (default value 2)

    FNMAI_SPARSE : as FNMAI plus
                'regul'  for l1-smoothness of X (default 0)

    NNSC     :  'alpha'       for stepsize of gradient update of A
                'sparse_par'  for sparsity

This module is based on:

    - Daniel D. Lee and H. Sebastian Seung:

          "Algorithms for non-negative matrix factorization", 
          in Advances in Neural Information Processing 13 
          (Proc. NIPS*2000) MIT Press, 2001.

          "Learning the parts of objects by non-negative matrix
           factorization",
          Nature, vol. 401, no. 6755, pp. 788-791, 1999.

    - A. Cichocki and A-H. Phan:

          "Fast local algorithms for large scale Nonnegative Matrix and
           Tensor Factorizations",
          IEICE Transaction on Fundamentals,
          in print March 2009.

    - P. O. Hoyer

          "Non-negative Matrix Factorization with sparseness
           constraints", 
          Journal of Machine Learning Research, vol. 5, pp. 1457-1469,
          2004.


    - Dongmin Kim, Suvrit Sra,Inderjit S. Dhillon:

           "Fast Newton-type Methods for the Least Squares Nonnegative Matrix
           Approximation Problem" 
           SIAM Data Mining (SDM), Apr. 2007 
          

    - Ngoc-Diep Ho:

        dissertation from
        http://edoc.bib.ucl.ac.be:81/ETD-db/collection/available/BelnUcetd-06052008-235205/
                

#############################################################################

Copyright (c) 2009 Uwe Schmitt, uschmitt@mineway.de

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
    * notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above
    * copyright notice, this list of conditions and the following
    * disclaimer in the documentation and/or other materials provided
    * with the distribution.  Neither the name of the <ORGANIZATION>
    * nor the names of its contributors may be used to endorse or
    * promote products derived from this software without specific
    * prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
import numpy
try:
    import mdp
    if mdp.__version__ >= '2.6':
        MDP = True
    else:
        MDP = False
except:
    MDP = False

from PyMca import py_nnma
DEBUG = 0

function_list = ['FNMAI', 'ALS', 'FastHALS', 'GDCLS']
function_dict = {"NNSC": py_nnma.NNSC,
                 "FNMAI_SPARSE": py_nnma.FNMAI_SPARSE,
                 "FNMAI": py_nnma.FNMAI,
                 "GDCLS_L1": py_nnma.GDCLS_L1,
                 "GDCLS": py_nnma.GDCLS,
                 "ALS": py_nnma.ALS,
                 "NMFKL": py_nnma.NMFKL,
                 "NMF": py_nnma.NMF,
                 "RRI": py_nnma.RRI,
                 "FastHALS": py_nnma.FastHALS,
                 "SNMF": py_nnma.SNMF,
                 }
def nnma(stack, ncomponents, binning=None,
         function=None, eps=5e-5, verbose=DEBUG, maxcount=1000, kmeans=False):
    if kmeans and (not MDP):
        raise ValueError("K Means not supported")
    #I take the defaults for the other parameters
    param = dict(alpha=.1, tau=2, regul=1e-2, sparse_par=1e-1, psi=1e-3)
    if function is None:
        function = 'FNMAI'
    nnma_function = function_dict[function]
    if binning is None:
        binning = 1

    if hasattr(stack, "info") and hasattr(stack, "data"):
        data = stack.data[:]
    else:
        data = stack[:]


    oldShape = data.shape

    if len(data.shape) == 3:
        r, c, N = data.shape
        if isinstance(data, numpy.ndarray):
            data.shape = r*c, N
    else:
        r, N = data.shape
        c = 1

    if isinstance(data, numpy.ndarray):
        if binning > 1:
            data=numpy.reshape(data,[data.shape[0], data.shape[1]/binning, binning])
            data=numpy.sum(data , axis=-1)
            N=N/binning
    else:
        oldData = data
        N = int(N/binning)
        try:
            data = numpy.zeros((r, c, N), oldData.dtype)
        except MemoryError:
            try:
                data = numpy.zeros((r, c, N), numpy.float32)
            except MemoryError:
                text  = "NNMAModule only works properly on numpy arrays.\n"
                text += "Memory Error: Higher binning may help."
                raise TypeError(text)
        if binning == 1:
            if len(oldShape) == 3:
                for i in range(r):
                    data[i,:,:] = oldData[i,:,:]
                data.shape = r * c, N
            else:
                data.shape = r * c, N
                for i in range(r*c):
                    data[i,:] = oldData[i,:]
        else:
            if len(oldShape) == 3:
                for i in range(r):
                    tmpData = oldData[i,:,:]
                    tmpData.shape = c, N, binning
                    data[i,:,:] = numpy.sum(tmpData, axis=-1)
                data.shape = r * c, N
            else:
                data.shape = r * c, N
                for i in range(r*c):
                    tmpData = oldData[i,:]
                    tmpData.shape = N, binning
                    data[i,:] =  numpy.sum(tmpData, axis=-1)

    #mindata = data.min()
    #numpy.add(data, -mindata+1, data)
    #I do not know the meaning of these paramenters
    #py_nnma.scale(newdata)
    param = dict(alpha=.1, tau=2, regul=1e-2, sparse_par=1e-1, psi=1e-3)
    #Start tolerance
    #1E+3 is conservative/fast
    #1E-3 is probably slow
    Astart = None
    Xstart = None
    #for i in range(start_ncomponents, ncomponents):
    converged = False
    while not converged:
        A, X, obj, count, converged = nnma_function(data,
                                                    ncomponents,
                                                    Astart,
                                                    Xstart,
                                                    eps=eps,
                                                    maxcount=maxcount,
                                                    verbose=verbose,
                                                    **param)
        if not converged:
            print("WARNING: Possible problems converging")
    #if binning > 1:
    #    numpy.add(data, mindata-1, data)
    #data.shape = oldShape
    images = A.T
    if 0:
        images.shape = ncomponents, r, c
        return images, numpy.ones((ncomponents), numpy.float32),X

    #order and scale images according to Gerd Wellenreuthers' recipe
    #normalize all maps to be in the range [0, 1]
    for i in range(ncomponents):
        norm_factor = numpy.max(images[i, :])
        if norm_factor > 0:
            images[i, :] *= 1.0/norm_factor
            X[i, :] *= norm_factor

    #sort NNMA-spectra and maps
    total_nnma_intensity = []
    for i in range(ncomponents):
        total_nnma_intensity += [[numpy.sum(images[i,:])*\
                                  numpy.sum(X[i,:]), i]]

    sorted_idx = [item[1] for item in sorted(total_nnma_intensity)]
    sorted_idx.reverse()

    #original data intensity
    original_intensity = numpy.sum(data)

    #final values
    if kmeans:
        n_more = 1
    else:
        n_more = 0
    new_images  = numpy.zeros((ncomponents + n_more, r*c), numpy.float32)
    new_vectors = numpy.zeros((X.shape[0]+n_more, X.shape[1]), numpy.float32)
    values      = numpy.zeros((ncomponents+n_more,), numpy.float32)
    for i in range(ncomponents):
        idx = sorted_idx[i]
        if 1:
            new_images[i, :] = images[idx, :]
        else:
            #imaging the projected sum gives same results
            Atmp = images[idx, :]
            Atmp.shape = -r*c, 1
            Xtmp = X[idx,:]
            Xtmp.shape = 1, -1
            new_images[i, :] = numpy.sum(numpy.dot(Atmp, Xtmp), axis=1)
        new_vectors[i,:] = X[idx,:]
        values[i] = 100.*total_nnma_intensity[idx][0]/original_intensity
    new_images.shape = ncomponents + n_more, r, c
    if kmeans:
        classifier = mdp.nodes.KMeansClassifier(ncomponents)
        for i in range(ncomponents):
            classifier.train(new_vectors[i:i+1])
        k = 0
        for i in range(r):
            for j in range(c):
                spectrum = data[k:k+1,:]
                new_images[-1, i,j] = classifier.label(spectrum)[0]
                k += 1
    return new_images, values, new_vectors        

if __name__ == "__main__":
    from PyMca import EDFStack
    from PyMca import EdfFile
    import os
    import sys
    import time
    inputfile = "D:\DATA\COTTE\ch09\ch09__mca_0005_0000_0000.edf"    
    if len(sys.argv) > 1:
        inputfile = sys.argv[1]
        print(inputfile)
    elif os.path.exists(inputfile):
        print("Using a default test case")
    else:
        print("Usage:")
        print("python NNMAModule.py indexed_edf_stack")
        sys.exit(0)
    stack = EDFStack.EDFStack(inputfile)
    r0, c0, n0 = stack.data.shape
    ncomponents = 10
    outfile = os.path.basename(inputfile)+"ICA.edf"
    e0 = time.time()
    images, eigenvalues, eigenvectors =  nnma(stack.data, ncomponents,
                                                     binning=1)
    print("elapsed = %f" % (time.time() - e0))
    if os.path.exists(outfile):
        os.remove(outfile)
    f = EdfFile.EdfFile(outfile)
    for i in range(ncomponents):
        f.WriteImage({}, images[i,:])
    sys.exit(0)
