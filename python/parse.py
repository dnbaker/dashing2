import scipy.sparse as sp
import numpy as np


def parse_knn(path, idsize=4, dstsize=4):
    '''

    Input: path [str]
           idsize=4
               Size of integers. Defaults to 4, can be changed with -DLSHIDTYPE=TYPE
               For instance, if compiled with -DLSHIDTYPE="uint64_t", you would need to use 8
           dstsize=4 
               Size of distance floating-point values. Defaults to 4 (floats), can be changed with -DDASHING2_INDEX_FLOAT_TYPE=TYPE
               For instance, if compiled with -DDASHING2_INDEX_FLOAT_TYPE="double", you would need to use 8
    Output: scipy.sparse.csr_matrix


    Dashing2's sparse notation is emitted in CSR (compressed-sparse row) notation.

    Format on disk:
    16 bytes: uint64_t nids, uint64_t nnz
    (nids + 1) * 8 bytes: indptr in uint64_t
    nnz * sizeof(LSHIDType): indices in LSHIDType (default uint32_t)
    nnz * sizeof(LSHDistType): data in LSHDistType (default float)
    '''
    import scipy.sparse as sp
    fp = open(path, "rb")
    ft = np.float32 if dstsize == 4 else np.float64
    if ft is np.float64:
        dstsize = 8
    it = {1: np.uint8, 2: np.uint16, 4: np.uint32,
          8: np.uint64}[idsize]
    nids, nnz = np.frombuffer(fp.read(16), dtype=np.uint64)
    indptr = np.frombuffer(fp.read(8 * (nids + 1)), dtype=np.uint64)
    indices = np.frombuffer(fp.read(idsize * nnz), dtype=it)
    data = np.frombuffer(fp.read(dstsize * nnz), dtype=ft)
    return sp.csr_matrix((data, indices, indptr), shape=(nids, nids))

__all__ = ["parse_knn"]
