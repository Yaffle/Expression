/*jshint esversion:6, bitwise:false*/

const BitSetWordSize = 31; // see https://v8.dev/blog/pointer-compression

function packedArray(n) {
  // `%DebugPrint(array)` in `node --allow-native-syntax`
  // see https://v8.dev/blog/elements-kinds
  const array = [];
  for (let i = 0; i < n; i++) {
    array.push(0);
  }
  return array.slice(0); // slice to reduce the size of the internal storage
}
function BitSet(size) {
  const n = Math.ceil(size / (4 * BitSetWordSize)) * 4;
  this.data = packedArray(n);
  this.size = size;
}
BitSet.prototype.nextSetBit = function (index) {
  if ((index | 0) >= (this.size | 0)) {
    return -1;
  }
  const data = this.data;
  let q = Math.floor(index / BitSetWordSize);
  let r = index % BitSetWordSize;
  let x = data[q] >> r;
  while (x === 0) {
    q += 1;
    if (q === data.length) {
      return -1;
    }
    x = data[q];
    r = 0;
  }
  if (x === (-1 << (BitSetWordSize - 1))) {
    // -x overflows
    r += BitSetWordSize - 1;
  } else {
    // https://stackoverflow.com/questions/61442006/whats-the-most-efficient-way-of-getting-position-of-least-significant-bit-of-a
    r += 31 - Math.clz32(x & -(+x));
  }
  return q * BitSetWordSize + r;
};
BitSet.prototype.toggle = function (index) {
  if ((index | 0) >= (this.size | 0)) {
    throw new RangeError();
  }
  const q = Math.floor(index / BitSetWordSize);
  const r = index % BitSetWordSize;
  this.data[q] ^= (r === BitSetWordSize - 1 ? ((-1) << r) : (1 << r));
};
BitSet.prototype.xor = function (other) {
  const a = this.data;
  const b = other.data;
  const n = a.length | 0;
  if (n !== b.length || n % 4 !== 0) {
    throw new RangeError();
  }
  for (let i = 0; i < n; i += 4) {
    a[i + 0] ^= b[i + 0] | 0;
    a[i + 1] ^= b[i + 1] | 0;
    a[i + 2] ^= b[i + 2] | 0;
    a[i + 3] ^= b[i + 3] | 0;
  }
};
BitSet.prototype.clear = function () {
  for (let i = 0; i < this.data.length; i += 1) {
    this.data[i] = 0;
  }
};
BitSet.prototype.toString = function () {
  return this.data.map(x => (x >>> 0).toString(2).padStart(BitSetWordSize, '0').split('').reverse().join('')).join('').slice(0, this.size);
};



// pass factorizations with associated values to the next call
// returns linear combinations of vectors which result in zero vector by modulo 2
// (basis of the kernel of the matrix)
function solve(matrixSize, bitsetRows = false) {
  // We build the augmented matrix in row-echelon form with permuted rows, which can grow up to matrixSize rows:
  // The augmented matrix is stored in the lower triangle!
  const M = new Array(matrixSize).fill(null); // We will fill the matrix so pivot elements will be placed on the diagonal
  const associatedValues = new Array(matrixSize).fill(undefined);
  let nextSolution = null;
  let state = 1;
  const iterator = {
    next: function solve(rawRowAndValue) {
      while (true) {
        if (state === 1) {
          state = 0;
          return {value: nextSolution, done: false};
        }
        state = 1;
        const rawRow = rawRowAndValue[0];
        const associatedValue = rawRowAndValue[1];
        let row = null;
        if (bitsetRows) {
          row = rawRow;
        } else {
          row = new BitSet(matrixSize);
          const reverseColumns = true; // makes it much faster when the data is more dense from the beginning (?)
          for (let j = 0; j < rawRow.length; j++) {
            const unitaryColumn = rawRow[j];
            const c = reverseColumns ? matrixSize - 1 - unitaryColumn : unitaryColumn;
            row.toggle(c);
          }
        }
        // add row to the matrix maintaining it to be in row-echelon form:
        for (let pivotColumn = row.nextSetBit(0); pivotColumn !== -1 && row != null; pivotColumn = row == null ? -1 : row.nextSetBit(pivotColumn + 1)) {
          const pivotRow = M[pivotColumn];
          if (pivotRow != null) {
            // row-reduction:
            row.xor(pivotRow);
            //console.assert(row.nextSetBit(pivotColumn) > pivotColumn || row.nextSetBit(pivotColumn) === -1);
            row.toggle(pivotColumn);
          } else {
            //row.toggle(matrixSize + pivotColumn);
            associatedValues[pivotColumn] = associatedValue;
            M[pivotColumn] = row;
            row = null;
          }
        }
        if (row != null) {
          // row has a solution
          // extract solution from the augmented part of the matrix:
          const solution = [];
          for (let i = row.nextSetBit(0); i !== -1; i = row.nextSetBit(i + 1)) {
            solution.push(associatedValues[i]);
          }
          solution.push(associatedValue);
          nextSolution = solution;
        } else {
          nextSolution = null;
        }
      }
      //console.log(M.filter(x => x != null).map(x => x.toString()).join('\n').replaceAll('0', ' '))
    }
  };
  iterator[globalThis.Symbol.iterator] = function () {
    return this;
  };
  return iterator;
}

// sparseColumns,
// build a matrix, do eliminations, ...
// keep indexes

function SparseGF2Matrix(n, m) {
  this.rows = new Array(n).fill(null);
  this.columns = new Array(m).fill(null).map(x => []);
}
SparseGF2Matrix.prototype.addRow = function (rawRow, rowIndex) {

  function removePairs(rawRow) {
    const row = rawRow.slice(0);
    row.sort((a, b) => a - b);
    let k = -1;
    for (let i = 0; i < row.length; i++) {
      if (i === 0 || row[k] !== row[i]) {
        k += 1;
        row[k] = row[i];
      } else {
        k -= 1;
      }
    }
    return row.slice(0, k + 1);
  }
  const row = removePairs(rawRow);
  
  this.rows[rowIndex] = row;
  for (let i = 0; i < row.length; i++) {
    this.columns[row[i]].push(rowIndex);
  }
};
SparseGF2Matrix.prototype.deleteColumn = function (columnIndex) {
  const column = this.columns[columnIndex];
  if (column != null) {
    for (let i = 0; i < column.length; i++) {
      let rowIndex = column[i];
      const k = this.rows[rowIndex].indexOf(columnIndex);
      console.assert(k !== -1);
      this.rows[rowIndex].splice(k, 1);
    }
    this.columns[columnIndex] = null;
  }
};

SparseGF2Matrix.prototype.deleteRow = function (rowIndex) {
  const row = this.rows[rowIndex];
  if (row != null) {
    for (let i = 0; i < row.length; i++) {
      const columnIndex = row[i];
      this.columns[columnIndex] = this.columns[columnIndex].filter(r => r != rowIndex);
    }
    this.rows[rowIndex] = null;
  }
};

SparseGF2Matrix.prototype.reduce = function (rowIndex, xors) {
  const row = this.rows[rowIndex];
  console.assert(row.length === 1);
  const e = row[0];
  const column = this.columns[e];
  for (let i = column.length - 1; i >= 0; i--) {
    if (column[i] !== rowIndex) {
      const j = column[i];
      const a = this.rows[j];
      const k = a.indexOf(e);
      console.assert(k !== -1);
      a.splice(k, 1);
      column.splice(i, 1);
      xors.push({a: j, b: rowIndex});
    }
  }
  console.assert(column.length === 1 && column[0] === rowIndex);
  this.columns[e] = null;
  this.rows[rowIndex] = null;
};


// structured Gaussian elimination:
// see The Factorization of the Ninth Fermat Number (Author(s): A. K. Lenstra, H. W. Lenstra, Jr., M. S. Manasse, J. M. Pollard), page 344
// and Solving Large Sparse Linear Systems Over Finite Fields (B. A. LaMacchia • A.M. Odlyzko)
function sparseSolve(columnsCount) {
  if (columnsCount < 1024 * 2) {
    return solve(columnsCount);
  }

  function getKthLargest(array, k) {
    array.sort((a, b) => a - b);
    return array[Math.max(array.length - k, 0)];
  }

  function deleteNonSparseColumns(sparseMatrix, count) {
    const counts = sparseMatrix.columns.map(c => c == null ? 0 : c.length);
    const t = getKthLargest(counts, count);
    for (let i = 0; i < sparseMatrix.columns.length; i++) {
      const column = sparseMatrix.columns[i];
      if (column != null && column.length >= t) {
        sparseMatrix.deleteColumn(i);
      }
    }
  }

  function eliminateRowsWithSingleColumn(sparseMatrix, eliminatedColumns) {
    const rows = [];
    for (let i = 0; i < sparseMatrix.rows.length; i++) {
      const row = sparseMatrix.rows[i];
      if (row != null && row.length === 1) {
        rows.push(i);
      }
    }
    for (let i = 0; i < rows.length; i++) {
      const rowIndex = rows[i];
      const row = sparseMatrix.rows[rowIndex];
      if (row != null && row.length === 1) {
        sparseMatrix.reduce(rowIndex, xors);
        eliminatedColumns[row[0]] = true;
      } else {
        console.assert(row != null && row.length === 0);
      }
    }
    return rows.length;
  }

  const M = [];
  const data = [];
  const xors = [];
  let newMatrix = null;
  let solutions = null;
  let lastJ = -1;

  function calculateNewMatrix() {
    const FIRST_SPARSE_PRIME_INDEX = 2585;

    const sparseMatrix = new SparseGF2Matrix(M.length, columnsCount);
    for (let i = 0; i < M.length; i++) {
      sparseMatrix.addRow(M[i].filter(p => p >= FIRST_SPARSE_PRIME_INDEX), i);
    }

    const eliminatedColumns = new Array(columnsCount).fill(false);

    const S = 128;
    let sparseColumns = -1;
    do {
      sparseColumns = sparseMatrix.columns.reduce((p, column) => p + (column == null ? 0 : 1), 0);
      //console.log('!sparseColumns:', sparseColumns);
      const c = Math.ceil(sparseMatrix.rows.length / S) - eliminateRowsWithSingleColumn(sparseMatrix, eliminatedColumns);
      if (c > 0) {
        deleteNonSparseColumns(sparseMatrix, c);
      }
    } while (sparseColumns > 0);

    const newColumnIds = new Array(columnsCount).fill(0);
    let k = 0;
    for (let i = 0; i < sparseMatrix.columns.length; i++) {
      if (!eliminatedColumns[i]) {
        newColumnIds[i] = k;
        k += 1;
      } else {
        newColumnIds[i] = -1;
      }
    }
    const newColumnsCount = k;

    const pool = [];
    function getNewRow(row) {
      let newRow = null;
      if (pool.length > 0) {
        newRow = pool.pop();
        newRow.clear();
      } else {
        newRow = new BitSet(newColumnsCount);
      }
      for (let j = 0; j < row.length; j++) {
        const e = newColumnIds[row[j]];
        if (e !== -1) {
          newRow.toggle(newColumnsCount - 1 - e);
        }
      }
      return newRow;
    }
    //newMatrix = M.map(row => getNewRow(row));

    //trying to reduce the memory usage by new rows:
    newMatrix = M.map(row => null);
    const lastUsage = sparseMatrix.rows.map(row => (row != null ? 1/0 : -1));
    for (let i = xors.length - 1; i >= 0; i -= 1) {
      if (lastUsage[xors[i].a] !== -1) {
        if (lastUsage[xors[i].b] === -1) {
          lastUsage[xors[i].b] = i;
        }
      } else {
        xors[i] = null;
      }
    }
    const newRow = function (i, j) {
      if (newMatrix[i] == null) {
        console.assert(lastUsage[i] >= j);
        newMatrix[i] = getNewRow(M[i]);
      }
      const result = newMatrix[i];
      return result;
    };

    for (let i = 0; i < xors.length; i++) {
      if (xors[i] != null) {
        newRow(xors[i].a, i).xor(newRow(xors[i].b, i));
        if (lastUsage[xors[i].b] <= i) {
          pool.push(newMatrix[xors[i].b]);
          newMatrix[xors[i].b] = null;
        }
      }
    }

    for (let i = 0; i < sparseMatrix.rows.length; i++) {
      if (sparseMatrix.rows[i] != null) {
        console.assert(sparseMatrix.rows[i].length === 0);
        newMatrix[i] = newRow(i, 1/0);
      } else {
        newMatrix[i] = null;
      }
    }

    /*console.time('n1');
    const template = packedArray(newColumnsCount);
    newMatrix = newMatrix.map(function (row) {
      if (row == null) {
        return null;
      }
      let k = 0;
      let from = -1;
      from = row.nextSetBit(from + 1);
      while (from !== -1) {
        template[k] = from;
        k += 1;
        from = row.nextSetBit(from + 1);
      }
      return template.slice(0, k);
    });
    console.timeEnd('n1');*/

    //console.log('newMatrix', newMatrix.filter(x => x != null).length);

    //const start = Date.now();
    solutions = solve(newColumnsCount, true); // find products of Y_k = Y, so that Y is a perfect square
    solutions.next();
  }

  let firstTime = true;
  const iterator = {
    next: function sparseSolve(rawRowAndValue) {
      if (firstTime) {
        firstTime = false;
        return {value: null, done: false};
      }
      if (solutions == null) {
        if (M.length < columnsCount + 32) {
          M.push(rawRowAndValue[0]);
          data.push(rawRowAndValue[1]);
          return {value: null, done: false};
        }
        calculateNewMatrix();
      }
      for (let j = lastJ + 1; j < newMatrix.length; j++) {
        lastJ = j;
        const row = newMatrix[j];
        if (row != null) {
          const solution = solutions.next([row, j]).value;
          if (solution != null) {
            //console.log('c', c.length);
            const set = new Array(M.length).fill(0);
            for (let i = 0; i < solution.length; i++) {
              if (set[solution[i]] !== 0) {
                throw new RangeError();
              }
              set[solution[i]] = 1;
            }
            for (let i = xors.length - 1; i >= 0; i--) {
              if (xors[i] != null && set[xors[i].a] === 1) {
                set[xors[i].b] = 1 - set[xors[i].b];
              }
            }
            //console.log('solution', solution, Array.from(s));
            const newSolution = [];
            for (let i = 0; i < set.length; i++) {
              if (set[i] === 1) {
                newSolution.push(data[i]);
              }
            }
            //console.log('a', Date.now() - start);
            return {value: newSolution, done: false};
          }
        }
      }
    }
  };
  iterator[globalThis.Symbol.iterator] = function () {
    return this;
  };
  return iterator;
}

solve.solve = solve;
solve.sparseSolve = sparseSolve;
sparseSolve.solve = solve;
sparseSolve.sparseSolve = sparseSolve;

export default solve;
