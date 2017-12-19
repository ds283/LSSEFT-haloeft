import numpy as np

import settings

# import a matrix from a flattened table
def import_matrix(matrix, table, tag, name):

    # generate flags matrix to catch import errors leading to uninitialized elements
    flags = np.zeros_like(matrix, dtype=bool)

    # populate matrix
    for row in table:

        i = row['i']-1
        j = row['j']-1

        matrix[i, j] = row['value']
        flags[i, j] = True

    # raise exception if any elements left uninitialized
    if np.any(flags == False):

        print '!! region {tag}: not all elements of {name} initialized'.format(tag=tag, name=name)
        print table
        raise RuntimeError


# invert a covariance matrix
def invert_covariance_matrix(matrix, tag, name):

    # get eigenvalues of matrix, needed to estimate condition number
    w, p = np.linalg.eig(matrix)

    # condition number is roughly ratio of largest to smallest eigenvalues
    condition_number = abs(np.min(w)/np.max(w))

    # if any negative eigenvalues or condition number is very small, use Moore-Penrose inverse
    if not np.all(w > 0) or condition_number < settings.ConditionNumberTolerance:

        print '-- region {tag}: using Moore-Penrose inverse for covariance matrix in "{name}" group'.format(tag=tag, name=name)

        inv = np.linalg.pinv(matrix)

    else:

        inv = np.linalg.inv(matrix)

    # validate that inverse matrix is numerically good
    test1 = abs(np.dot(matrix, inv) - np.identity(matrix.shape[0]))
    test2 = abs(np.dot(inv, matrix) - np.identity(matrix.shape[0]))

    if test1.max() > settings.MatrixInverseTolerance or test2.max() > settings.MatrixInverseTolerance:
        print '-- region {tag}: test1 max value = {t1}, test2 max value = {t2}'.format(tag=tag, t1=test1.max(), t2=test2.max())
        raise RuntimeError

    return inv
