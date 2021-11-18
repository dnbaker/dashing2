import scipy.sparse as sp
import numpy as np
from collections import namedtuple

ParsedSignatureMatrix = namedtuple("ParsedSignatureMatrix", 'nseqs, cardinalities, signatures')
ParsedKmerMatrix = namedtuple("ParsedKmerMatrix", 'k,w,canon,alphabet,sketchsize,seed,kmers')


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



def parse_binary_signatures(path):
    '''
    If -o is specified for dashing2, this file has signatures and metadata saved to disk.
    Parsing this with `parse_binary_signatures` will yield a ParsedSignatureMatrix tuple
    consisting of (nseqs, cardinalities, signatures).
    '''
    dat = np.memmap(path, np.uint8) # Map whole file
    nseqs, sketchsize = map(int, dat[:16].view(np.uint64))
    cardinalities = dat[16:16 + (8 * nseqs)].view(np.float64)
    signatures = dat[16 + (8 * nseqs):].view(np.float64).reshape(nseqs, -1)
    sigmul = sketchsize // signatures.shape[1]
    if sigmul != 1:
        signatures = signatures.view({2: np.uint32, 1: np.uint64, 4: np.uint16, 8: np.uint8}[sigmul])
    return ParsedSignatureMatrix(nseqs, cardinalities, signatures)


def parse_binary_kmers(path):
    '''
    If --save-kmers is specified for dashing2, this file has kmers and metadata saved to disk.
    Parsing this with `parse_binary_kmers` will yield a ParsedKmerMatrix tuple
    consisting of (k, w, canon, alphabet, sketchsize, signatures).
    (nseqs is implicit in the matrix because sketchsize is specified.)
    '''
    dat = np.memmap(path, np.uint8)
    d, s, k, w  = map(int, dat[:16].view(np.uint32))
    canon = bool((d >> 8) & 1)
    alph = alphabetcvt[d & 0xff]
    seed = int(dat[16:24].view(np.uint64))
    kmers = dat[24:].view(np.uint64).reshape(-1, s)
    return ParsedKmerMatrix(k, w, canon, alph, s, seed, kmers)


def alphabetcvt(x):
    itd = {"DNA": 0, "BYTES": 1, "PROTEIN20": 2, "PROTEIN14": 3, "PROTEIN6": 4, "DNA2": 5, "DNAC": 6}
    itl = ["DNA", "BYTES", "PROTEIN20", "PROTEIN14", "PROTEIN6", "DNA2", "DNAC"]
    if isinstance(x, str):
        return itd[x.upper()]
    if isinstance(x, int):
        return itl[x]
    raise InvalidArgument("alphabetcvt only accepts strings (which are decoded to integers to match Bonsai's alphabets) or integers (which returns a string describing the alphabet.")


def pairwise_equality_compare(input_matrix, nthreads=1):
    '''
    This performs pairwise equality comparisons along a signature matrix.
    Inputs:
        input_matrix - 2D Numpy Matrix of dimensions (nrec, nsigs)
        nthreads = 1: number of threads to use. This is only used if sketch (https://github.com/dnbaker/sketch) has been installed.

    If `sketch` is not installed, computation will be performed serially via NumPy rather than using parallelized and SIMD-accelerated comparisons.
    Return a distance matrix of size nrec-choose-2 (nrec * (nrec - 1) // 2).
    This can be converted into a complete distance matrix using scipy.spatial.distance.squareform.
    '''
    assert isinstance(input_matrix, np.ndarray), "Expected a numpy array for this function."
    assert input_matrix.ndim == 2, "Expected a 2d array"
    nthreads = max(nthreads, 1)
    try:
        import sketch
        return sketch.pcmp(input_matrix, nthreads=nthreads)
    except ImportError:
        print("sketch python library not installed (see https://github.com/dnbaker/sketch). Falling back to numpy computation.", file=sys.stderr)
        nr, nc = input_matrix.shape
        dt = np.uint8 if nr <= 0xFF else (np.uint16 if nr <= 0xFFFF else (np.uint32 if nr <= 0xFFFFFFFF else np.uint64))
        shape = nr * (nr - 1) // 2
        ret = np.zeros(shape, dtype=dt)
        idx = 0
        for i in range(nr):
            lc = nr - i - 1
            ret[idx:idx + lc] = np.sum(input_matrix[i] == input_matrix[i + 1:nr], axis=1)
            idx += lc
        return ret


def parse_binary_clustering(path, d64=False):
    '''
    Parses binary output for dashing2 --greedy computation.
    Yields a list of numpy arrays of integers specifying which entities belong in which clusters, one for each cluster.
    These clusters are fully distinct.
    '''
    data = np.memmap(fpath, np.uint8)
    nclusters, nsets = map(int, data[:16].view(np.uint64))
    endp = 16 + 8 * nclusters
    indptr = data[16:endp].view(np.uint64)
    indices = data[endp:].view(np.uint64 if d64 else np.uint32)
    return [indices[start:end] for start, end in zip(indptr[:-1], indptr[1:])]


def parse_binary_distmat(path):
    '''
    Parse all-pairs distances from binary distance matrix at <path>.
    '''
    return np.memmap(path, np.float32)


def parse_binary_rectmat(path, fpath, qpath):
    '''
    Parse distance matrix from path using <fpath> and <qpath> as reference/query pairs.
    fpath must have been provided to Dashing2 with -F/--ffile
    qpath must have been provided to Dashing2 with -Q/--qfile
    '''
    ffiles, qfiles = map(lambda x: list(map(str.strip, open(x))), (fpath, qpath))
    nref, nquery = map(len, (ffiles, qfiles))
    return np.memmap(path, np.float32).reshape(nref, nquery)


def parse_binary_contain(path):
    rawdata = np.memmap(path, dtype=np.float32)
    nref, nqueries = map(int, rawdata[:4].view(np.uint64))
    coverage_fractions = rawdata[4:4 + nref * nqueries].reshape(nqueries, nref)
    meandepth = rawdata[4 + nref * nqueries:].reshape(nqueries, nref)
    return {"nref": nref, "nqueries": nqueries, "coverage_matrix": coverage_fractions, "depth_matrix": meandepth}


__all__ = ["parse_knn", "parse_binary_signatures", "ParsedSignatureMatrix", "parse_binary_kmers", "ParsedKmerMatrix", "alphabetcvt", "pairwise_equality_compare", "parse_binary_clustering", "parse_binary_distmat", "parse_binary_rectmat"]
