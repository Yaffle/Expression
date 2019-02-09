/*global Matrix, Expression*/

(function () {
"use strict";

// https://ca.wikipedia.org/wiki/Forma_canònica_de_Jordan

// https://es.wikipedia.org/wiki/Forma_canónica_de_Jordan

Expression.getFormaDeJordan = function (matrix, eigenvalues, multiplicities, hack) {
  function getSolutionSet(matrix) {
    var fullMatrix = matrix.augment(Matrix.Zero(matrix.cols(), 1));
    //TODO: Matrix.GaussMontante
    var result = fullMatrix.toRowEchelon(Matrix.GaussJordan, "solving", undefined);
    var tmp = Matrix.solveByGaussNext(result.matrix);
    var currentEigenvectors = Matrix.getSolutionSet(tmp).basisVectors;
    return currentEigenvectors;//?
  }
  function matrixFromBlocks(blocks) {
    var start = 0;
    var J = Matrix.Zero(n, n);
    for (var i = 0; i < blocks.length; i += 1) {
      var b = blocks[i];
      J = J.map(function (e, i, j) {
        if (i >= start && i < start + b.size) {
          return i === j ? b.eigenvalue : (i !== start + b.size - 1 && j === i + 1 ? Expression.ONE : Expression.ZERO);
        }
        return e;
      });
      start += b.size;
    }
    return J;
  }
  function matrixFromBasis(basis) {
    if (basis.length === 0) {
      throw new Error();
    }
    return Matrix.Zero(basis[0].rows(), basis[0].rows()).map(function (e, i, j) {
      return basis[i].e(j, 0);
    });
  }
  function isSolution(coefficientMatrix, vector) {
    var f = coefficientMatrix.multiply(vector);
    if (f.cols() !== 1) {
      throw new RangeError("assertion failed");
    }
    for (var i = 0; i < f.rows(); i += 1) {
      if (!f.e(i, 0).equals(Expression.ZERO)) {
        return false;
      }
    }
    return true;
  }

  //!TODO: remove
  if (eigenvalues.length === matrix.rows()) {
    var eigenvectors = Expression.getEigenvectors(matrix, eigenvalues).eigenvectors;
    var tmp = Expression.diagonalize(matrix, eigenvalues, multiplicities, eigenvectors);
    var P = tmp.T;
    var J = tmp.L;
    var P_INVERSED = tmp.T_INVERSED;
    //console.log("P=" + P.toString() + ", J=" + J.toString());
    return {
      P: P,
      J: J,
      P_INVERSED: P_INVERSED
    };
  }
  //!

  var A = matrix;
  var n = A.rows();

  var basis = [];
  var blocks = [];
  //var bs = null;
  for (var i = 0; i < eigenvalues.length; i += 1) {
    // https://en.wikipedia.org/wiki/Generalized_eigenvector#Computation_of_generalized_eigenvectors
    var eigenvalue = eigenvalues[i];
    var algebraicMultiplicity = multiplicities[i];
    var B = A.subtract(Matrix.I(n).scale(eigenvalue));
    var m = 1;
    while (B.pow(m).rank() > n - algebraicMultiplicity) {
      m += 1;
    }
    m += 1;
    while (--m >= 1) {
      var z = 0;
      var pm = B.pow(m - 1).rank() - 2 * B.pow(m).rank() + B.pow(m + 1).rank();
      var solutionSet = getSolutionSet(B.pow(m));  // "kernel of A"
      for (var j = 0; j < solutionSet.length && z < pm; j += 1) {
        var solution = solutionSet[j];
        //console.log(B.pow(m).augment(solution).rank(), m, n);
        // (bs == null || bs.augment(solution).rank() > bs.rank())
        if (!isSolution(B.pow(m - 1), solution)) {
          var tmp = [];
          var s = solution;
          for (var k = 0; k < m; k += 1) {
            tmp.push(s);
            //bs = bs == null ? s : bs.augment(s);
            s = B.multiply(s);
          }
          z += 1;
          tmp.reverse();
          basis = basis.concat(tmp);
          blocks.push({
            size: m,
            eigenvalue: eigenvalue
          });
        }
      }
    }
  }
  var J = matrixFromBlocks(blocks);
  if (basis.length !== n) {
    throw new RangeError("assertion failed");
  }
  var P = matrixFromBasis(basis).transpose();
  //console.log("P=" + P.toString() + ", J=" + J.toString());
  //var P_INVERSED = P.inverse();
  var P_INVERSED = P.isExact() ? P.inverse() : (hack ? null : getInverse(A, eigenvalues, multiplicities, P));
  if (!hack && P.isExact()) {
  if (A.toString() !== P.multiply(J).multiply(P_INVERSED).toString()) {
    throw new RangeError("assertion failed");
  }
  }
  return {
    P: P,
    J: J,
    P_INVERSED: P_INVERSED
  };
};

// A = P*J*P^-1
// A^T = (P^-1)^T*J^T*P^T
// Note:
// (0 0 0 1)         (0 0 0 1)
// (0 0 1 0)         (0 0 1 0)
// (0 1 0 0) * J^T * (0 1 0 0) = J
// (1 0 0 0)         (1 0 0 0)
// where on the left we are doing row spaws, then doing column swaps.
// Note: the inverse of the anti-diagonal unitary matrix is the matrix itself.
// A^T = X*Y*X^-1
// (P^-1)^T*J^T*P^T = X*Y*X^-1
// (P^-1)^T*B^-1*J*B*P^T = X*Y*X^-1
// Then P^-1 = (X*B)^T .

var getInverse = function (A, eigenvalues, multiplicities, P) {
  // https://en.wikipedia.org/wiki/Diagonalizable_matrix : The row vectors of P^−1 are the left eigenvectors of A
  // https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Left_and_right_eigenvectors :  a left eigenvector of A is the same as the transpose of a right eigenvector of A^T, with the same eigenvalue
  var AT = A.transpose();
  var tmp2 = Expression.getFormaDeJordan(AT, eigenvalues, multiplicities, true);
  var J = tmp2.J;
  var X = tmp2.P;

  var n = A.cols();
  var B = Matrix.Zero(n, n).map(function (e, i, j) {
    function getCurrentBlock() {
      var s = i;
      while (s - 1 >= 0 && s < n && J.e(s - 1, s).equals(Expression.ONE)) {
        s -= 1;
      }
      var e = i;
      while (e + 1 < n && J.e(e, e + 1).equals(Expression.ONE)) {
        e += 1;
      }
      return {s: s, e: e};
    }
    var tmp = getCurrentBlock();
    return tmp.s + tmp.e === i + j ? Expression.ONE : Expression.ZERO;
  });

  var P_INVESRED = X.multiply(B).transpose();
  return Expression._unscaleInverseMatrix(P_INVESRED, P);
};

}());

