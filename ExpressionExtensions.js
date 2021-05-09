import Expression from './Expression.js';
import Matrix from './Matrix.js';
import Polynomial from './Polynomial.js';

Expression.getPolynomialRootsWithSteps = function (polynomial, fractionDigits, callback) {
  var roots = polynomial.getroots(callback);

  //TODO: tests
  //!2018-05-28
  //!2018-07-11
  // experimental code
  var zeros = [];
  if (typeof polynomial.getZeros === "function" && roots.length !== polynomial.getDegree()) {
    var p = Polynomial.of(Expression.ONE);
    for (var i = 0; i < roots.length; i += 1) {
      p = p.multiply(Polynomial.of(roots[i].negate(), Expression.ONE));
    }
    var r = polynomial.divideAndRemainder(p).quotient;
    var precision = Math.max(fractionDigits || 0, 5);
    zeros = r.getZeros(precision, true);
  }
  //!

  // removing of duplicates
  var uniqueRoots = [];
  var multiplicities = [];
  for (var i = 0; i < roots.length; i += 1) {
    var root = roots[i];
    var isDuplicate = false;
    var j = -1;
    while (++j < uniqueRoots.length) {
      if (uniqueRoots[j].equals(root)) {
        isDuplicate = true;
        multiplicities[j] += 1;
      }
    }
    if (!isDuplicate) {
      uniqueRoots.push(root);
      multiplicities.push(1);
    }
  }

  var m = 0;
  for (var i = 0; i < zeros.length; i += 1) {
    m += 1;
    var zero = zeros[i];
    var next = i + 1 < zeros.length ? zeros[i + 1] : undefined;
    if (next !== zero) {
      uniqueRoots.push(zero);
      multiplicities.push(m);
      m = 0;
    }
  }

  return {
    uniqueRoots: uniqueRoots,
    multiplicities: multiplicities
  };
};

Expression.getEigenvalues = function (matrix, fractionDigits, callback) {

  if (!matrix.isSquare()) {
    throw new RangeError("NonSquareMatrixException");
  }
  // TODO: remove Polynomial

  var determinant = matrix.map(function (e, i, j) {
    var p = i === j ? Polynomial.of(e, Expression.ONE.negate()) : (e.equals(Expression.ZERO) ? Polynomial.ZERO : Polynomial.of(e));
    return new Expression.Polynomial(p);
  }).determinant();
  determinant = determinant.polynomial;

  //!new (sin/cos)
  //TODO: fix
  determinant = determinant.map(function (e) { return e.simplifyExpression(); });

  var characteristicPolynomial = determinant;//!TODO: fix

  var tmp = Expression.getPolynomialRootsWithSteps(characteristicPolynomial, fractionDigits, callback);
  var uniqueRoots = tmp.uniqueRoots;
  var multiplicities = tmp.multiplicities;

  var eigenvalues = uniqueRoots;
  return {
    characteristicPolynomial: characteristicPolynomial,
    eigenvalues: eigenvalues,
    multiplicities: multiplicities
  };
};

Expression.getEigenvectors = function (matrix, eigenvalues) {
  var eigenvectors = [];
  for (var i = 0; i < eigenvalues.length; i += 1) {
    var n = matrix.cols();
    // matrix - E * eigenvalue
    var fullMatrix = matrix.subtract(Matrix.I(n).scale(eigenvalues[i])).augment(Matrix.Zero(n, 1));
    //TODO: Matrix.GaussMontante
    var result = fullMatrix.toRowEchelon(Matrix.GaussJordan, "solving", undefined);
    var tmp = Matrix.solveByGaussNext(result.matrix);
    var currentEigenvectors = Matrix.getSolutionSet(tmp).basisVectors;
    eigenvectors = eigenvectors.concat(currentEigenvectors);//?
  }
  return {
    eigenvectors: eigenvectors
  };
};

var getInverse = function (A, eigenvalues, T) {
  // https://en.wikipedia.org/wiki/Diagonalizable_matrix : The row vectors of P^−1 are the left eigenvectors of A
  // https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Left_and_right_eigenvectors :  a left eigenvector of A is the same as the transpose of a right eigenvector of A^T, with the same eigenvalue
  var AT = A.transpose();
  var tmp2 = Expression.getEigenvectors(AT, eigenvalues);
  var eigenvectors = tmp2.eigenvectors;
  var T_INVERSED = Matrix.I(AT.cols()).map(function (e, i, j) {
    return eigenvectors[i].e(j, 0);
  });
  return _unscaleInverseMatrix(T_INVERSED, T);
};

var _unscaleInverseMatrix = function (T_INVERSED, T) {
  // we know, that the result is {{s_1, 0, 0, 0}, {0, s_2, 0, 0}, {0, 0, s_3, 0}, {0, 0, 0, s_4}}
  var trickyMultiply = function (a, b) {
    var n = a.rows();
    return Matrix.Zero(n, n).map(function (element, i, j) {
      if (i !== j) {
        return Expression.ZERO;
      }
      var rows = n;
      var k = -1;
      while (++k < rows) {
        var current = a.e(i, k).multiply(b.e(k, j));
        element = k === 0 ? current : element.add(current);
      }
      return element;
    });
  };

  var S = trickyMultiply(T_INVERSED, T);
  var S_INVERSED = S.map(function (e, i, j) {
    return i === j ? e.inverse() : Expression.ZERO;
  });

  return S_INVERSED.multiply(T_INVERSED);
};
Expression._unscaleInverseMatrix = _unscaleInverseMatrix;//TODO: make private

// A = T^-1 L T ,T-matrix of own vectors, L - matrix of own values

Expression.diagonalize = function (matrix, eigenvalues, multiplicities, eigenvectors) {
  if (!matrix.isSquare()) {
    throw new RangeError("NonSquareMatrixException");
  }
  if (Expression.sum(multiplicities) !== matrix.cols()) {
    throw new RangeError();
  }
  if (eigenvectors.length !== matrix.cols()) {
    throw new RangeError();
  }
  // https://en.wikipedia.org/wiki/Jordan_normal_form
  // A is diagonalizable if and only if, for every eigenvalue λ of A, its geometric and algebraic multiplicities coincide.

  // TODO: text
  //!!!
  var getEigenvalue = function (i) {
    var eigenvalueIndex = -1;
    var s = i;
    while (s >= 0) {
      s -= multiplicities[eigenvalueIndex + 1];
      eigenvalueIndex += 1;
    }
    return eigenvalues[eigenvalueIndex];
  };
  var L = Matrix.I(matrix.cols()).map(function (element, i, j) {
    return (i === j ? getEigenvalue(i) : Expression.ZERO);
  });
  var T = Matrix.I(matrix.cols()).map(function (e, i, j) {
    return eigenvectors[j].e(i, 0);
  });

  //var T_INVERSED = T.inverse();
  var T_INVERSED = T.isExact() ? T.inverse() : getInverse(matrix, eigenvalues, T);

  return {T: T, L: L, T_INVERSED: T_INVERSED};
};

Expression.LUDecomposition = function (matrix) {
  //https://en.wikipedia.org/wiki/LU_decomposition#Code_examples
  //TODO: remove(?) - matrix.toRowEchelon(...)
  var N = matrix.rows();
  var a = matrix;
  var Lower = Matrix.I(N);
  var P = Matrix.I(N);
  var swapFlag = false;
  var pivotRow = 0;
  for (var n = 0; n < matrix.cols(); n += 1) {
    if (pivotRow < N) {
      var c = pivotRow;
      if (a.e(pivotRow, n).equals(Expression.ZERO)) {
        for (var k = pivotRow + 1; k < N && c === pivotRow; k += 1) {
          if (!a.e(k, n).equals(Expression.ZERO)) {
            c = k;
          }
        }
        if (c !== pivotRow) {
          var S = Matrix.I(N);
          S = S.map(function (element, i, j) {
            return i === pivotRow ? S.e(c, j) : (i === c ? S.e(pivotRow, j) : element);
          });
          a = S.multiply(a);
          Lower = S.multiply(Lower.subtract(Matrix.I(N))).add(Matrix.I(N));
          P = S.multiply(P);
          swapFlag = true;
        }
      }
      if (!a.e(pivotRow, n).equals(Expression.ZERO)) {
        var L = Matrix.I(N).map(function (element, i, j) {
          return j === pivotRow && i >= pivotRow + 1 ? a.e(i, n).divide(a.e(pivotRow, n)).negate() : element;
        });
        a = L.multiply(a);
        Lower = Lower.multiply(L);
        pivotRow += 1;
      }
    }
  }
  Lower = Lower.map(function (element, i, j) {
    return i === j ? element : element.negate();
  });
  return {
    swapFlag: swapFlag,
    P: new Expression.Matrix(P),
    A: new Expression.Matrix(matrix),
    L: new Expression.Matrix(Lower),
    U: new Expression.Matrix(a)
  };
};

Expression.isReal = function (e) {
  var isReal = function (e) {
    if (e instanceof Expression.Integer) {
      return true;
    }
    if (e instanceof Expression.NthRoot) {
      return isReal(e.a);
    }
    if (e instanceof Expression.BinaryOperation) {
      return isReal(e.a) && isReal(e.b);
    }
    if (e === Expression.E || e === Expression.PI) {
      return true;
    }
    if (e instanceof Expression.Function) {
      return isReal(e.a);
    }
    if (e instanceof Expression.PolynomialRootSymbol) {
      return true;//TODO: ?
    }
    if (e instanceof Expression.ExpressionWithPolynomialRoot) {
      return isReal(e.e);
    }
    if (e instanceof Expression.ExpressionPolynomialRoot) {
      return true;
    }
    return false;
  };
  return isReal(e);
};
Expression.isRealMatrix = function (A) {
  for (var i = 0; i < A.rows(); i += 1) {
    for (var j = 0; j < A.cols(); j += 1) {
      if (!Expression.isReal(A.e(i, j))) {
        return false;
      }
    }
  }
  return true;
};
Expression.CholeskyDecomposition = function (matrix) {
  var A = matrix;

  // check if A is square
  if (!A.isSquare()) {
    throw new RangeError("NonSquareMatrixException");
  }

  var n = A.rows();

  // check if A from R
  var isReal = Expression.isRealMatrix(A);

  // check if A is symmetric
  for (var i = 0; i < n; i += 1)  {
    for (var j = i; j < n; j += 1) {
      if (!A.e(i, j).equals(A.e(j, i).complexConjugate())) {
        if (isReal) {
          throw new RangeError("NonSymmetricMatrixException");
        } else {
          throw new RangeError("NonHermitianMatrixException");
        }
      }
    }
  }

  var L = new Array(n);
  for (var i = 0; i < n; i += 1) {
    L[i] = new Array(n);
    for (var j = 0; j < n; j += 1) {
      L[i][j] = Expression.ZERO;
    }
  }

  for (var j = 0; j < n; j += 1) {
    for (var i = j; i < n; i += 1) {
      var e = null;
      if (j === i) {
        var sum = null;
        for (var k = 0; k < j; k += 1) {
          var s = L[j][k].multiply(L[j][k].complexConjugate());
          sum = sum == null ? s : sum.add(s);
        }
        var x = sum == null ? A.e(j, j) : A.e(j, j).subtract(sum);
        //TODO: fix
        if (!Expression._isPositive(x)) {
          throw new RangeError("NonPositiveDefiniteMatrix");
        }
        e = x.squareRoot();
      } else {
        var sum = null;
        for (var k = 0; k < j; k += 1) {
          var x = L[i][k].multiply(L[j][k].complexConjugate());
          sum = sum == null ? x : sum.add(x);
        }
        e = (sum == null ? A.e(i, j) : A.e(i, j).subtract(sum)).divide(L[j][j]);
      }
      L[i][j] = e;
      console.log("l_%d%d = %s", i + 1, j + 1, L[i][j].toString());
    }
  }
  return {
    L: Matrix.padRows(L, null)
  };
};


Matrix.prototype.conjugateTranspose = function () {
  return this.transpose().map(e => e.complexConjugate());
};

// SVD-decomposition
Expression.SVDDecomposition = function (matrix) {
  // TODO: see email 
  // https://en.wikipedia.org/wiki/Singular_value_decomposition#Calculating_the_SVD
  // TODO: see https://web.mit.edu/be.400/www/SVD/Singular_Value_Decomposition.htm
  // The left-singular vectors of M are a set of orthonormal eigenvectors of MM*.
  var helper = function (matrix, eigenvalues, multiplicities) {
    var tmp = Expression.getEigenvectors(matrix, eigenvalues);
    var eigenvectors = tmp.eigenvectors;
    //TODO: orthogonalization - ?
    if (multiplicities.length !== matrix.rows()) {//TODO: https://math.stackexchange.com/questions/82467/eigenvectors-of-real-symmetric-matrices-are-orthogonal#answer-82471
      eigenvectors = GramSchmidtOrthogonalization(eigenvectors);
    }
    var norm = function (vector) {
      return vector.dot(vector).squareRoot();
    };
    //TODO: rewrite (use array of vectors - ?)
    var normInverses = new Array(eigenvectors.length).fill(null).map((e, j) => Expression.ONE.divide(norm(eigenvectors[j])));
    var eigenvectorMatrix = Matrix.Zero(eigenvectors[0].dimensions(), eigenvectors.length).map((e, i, j) => eigenvectors[j].e(i).multiply(normInverses[j]));
    return eigenvectorMatrix;
  };
  var MMstar = matrix.multiply(matrix.conjugateTranspose());
  var tmp = Expression.getEigenvalues(MMstar);
  var eigenvalues = tmp.eigenvalues;
  //!
  eigenvalues = eigenvalues.map(eigenvalue => eigenvalue instanceof Expression.ExpressionWithPolynomialRoot || eigenvalue instanceof Expression.ExpressionPolynomialRoot ? eigenvalue.upgrade() : eigenvalue);
  //!
  eigenvalues.sort((a, b) => a._pow(2).subtract(b._pow(2)).compareTo(Expression.ZERO) > 0 ? -1 : 1);
  var multiplicities = tmp.multiplicities;
  var U = helper(MMstar, eigenvalues, multiplicities);
  var diagonal = [];
  for (var i = 0; i < eigenvalues.length && diagonal.length < Math.min(matrix.rows(), matrix.cols()); i += 1) {
    var eigenvalue = eigenvalues[i];
    var entry = eigenvalue.squareRoot();
    var eigenvalue = eigenvalues[i];
    // https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Eigenspaces,_geometric_multiplicity,_and_the_eigenbasis_for_matrices
    //TODO: optimize - ?
    var geometricMultiplicity = multiplicities[i] === 1 ? 1 : MMstar.cols() - MMstar.subtract(Matrix.I(MMstar.cols()).scale(eigenvalue)).rank();//?
    for (var j = 0; j < geometricMultiplicity; j += 1) {
      diagonal.push(entry);
    }
  }
  while (diagonal.length < Math.min(matrix.rows(), matrix.cols())) {//?TODO: ?
    diagonal.push(Expression.ZERO);
  }
  var Sigma = Matrix.Zero(matrix.rows(), matrix.cols()).map((e, i, j) => (i === j ? diagonal[i] : Expression.ZERO));
  //var Sigma = Matrix.Diagonal(diagonal);
  //if (matrix.cols() > matrix.rows()) {//TODO: ?
  //  Sigma = Sigma.augment(Matrix.Zero(matrix.rows(), matrix.cols() - matrix.rows()));
  //}
  //var Sigma = Matrix.Zero(matrix.rows(), matrix.cols()).map((e, i, j) => (i === j ? eigenvalues[i].squareRoot() : Expression.ZERO));
  //var Vstar = ExpressionParser.parse(matrix.toString()).transformEquality(ExpressionParser.parse(U.multiply(Sigma).toString() + '*' + 'X', ExpressionParser.parse.c).simplify());
  
  //var s = Matrix.Zero(matrix.rows(), matrix.rows()).map((e, i, j) => (i === j ? (i < eigenvalues.length && !eigenvalues[i].equals(Expression.ZERO) ? eigenvalues[i].squareRoot().inverse() : Expression.ONE) : Expression.ZERO));
  var MstarM = matrix.conjugateTranspose().multiply(matrix);
  //var tmp = Expression.getEigenvalues(MstarM);
  //var eigenvalues = tmp.eigenvalues;
  //var multiplicities = tmp.multiplicities;
  var Vstar = helper(MstarM, eigenvalues, multiplicities).conjugateTranspose();
  //var Vstar = s.multiply(U.inverse().multiply(matrix));
  //var Vstar = U.conjugateTranspose().multiply(matrix);
  //console.log(U.multiply(U.conjugateTranspose()).toString());
  //console.log(Vstar.multiply(Vstar.conjugateTranspose()).toString());
  return {U: U, Sigma: Sigma, Vstar: Vstar};
};
