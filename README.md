
It is a homemade library for JavaScript.
It can parse expressions, solve and simplify systems of linear equations, find eigenvalues and eigenvectors.

Installation
============
npm install Yaffle/Expression

Usage example
=============
```html
<script src="./BigIntWrapper.js"></script>
<script src="./BigInteger.js"></script>
<script src="./BigIntegerMath.js"></script>
<script src="./Expression.js"></script>
<script src="./toDecimalString.js"></script>
<script src="./GF2.js"></script>
<script src="./sin.js"></script>
<script src="./complex.js"></script>
<script src="./NonSimplifiedExpression.js"></script>
<script src="./ExpressionParser.js"></script>
<script src="./Polynom.js"></script>
<script src="./Matrix.js"></script>
<script src="./toLaTeX.js"></script>
<script src="./polynomial-roots-finding.js"></script>
<script src="./FormaCanonicalDeJordan.js"></script>
<script src="./ExpressionExtensions.js"></script>
<script src="./Condition.js"></script>
<script src="./toMathML.js"></script>
<script>

  var p = Polynom.toPolynomial(ExpressionParser.parse("10x^5−17x^4−505x^3+1775x^2−249x−630"), ExpressionParser.parse("x"));
  console.log(p.getroots().toString()); // -1/2,5,21/5,(-73^0.5-7)/2,(73^0.5-7)/2

  var p = Polynom.toPolynomial(ExpressionParser.parse("x^5−2x^4−11x^3+26x^2−2x−13"), ExpressionParser.parse("x"));
  console.log(p.getZeros(5).result.toString()); // -3.412,-0.609,1.075,1.925,3.021

  var matrix = ExpressionParser.parse('{{1,2,3},{4,5,6},{7,8,9}}').matrix;
  console.log('matrix: ' + matrix.toString()); // matrix: {{1,2,3},{4,5,6},{7,8,9}}

  var x = Expression.getEigenvalues(matrix);
  var multiplicities = x.multiplicities;
  var eigenvalues = x.eigenvalues;
  console.log('eigenvalues: ' + eigenvalues.toString()); // eigenvalues: 0,(-3*33^0.5+15)/2,(3*33^0.5+15)/2

  var eigenvectors = Expression.getEigenvectors(matrix, x.eigenvalues).eigenvectors;
  console.log('eigenvectors: ' + eigenvectors.toString()); // eigenvectors: {{1},{-2},{1}},{{(-3*33^0.5-11)/22},{(-3*33^0.5+11)/44},{1}},{{(3*33^0.5-11)/22},{(3*33^0.5+11)/44},{1}}

  var y = Expression.diagonalize(matrix, eigenvalues, multiplicities, eigenvectors);
  console.log('diagonalization: ' + matrix.toString() + ' = ' + y.T.toString() + " * " + y.L.toString() + " * " + y.T_INVERSED.toString()); // diagonalization: {{1,2,3},{4,5,6},{7,8,9}} = {{1,(-3*33^0.5-11)/22,(3*33^0.5-11)/22},{-2,(-3*33^0.5+11)/44,(3*33^0.5+11)/44},{1,1,1}} * {{0,0,0},{0,(-3*33^0.5+15)/2,0},{0,0,(3*33^0.5+15)/2}} * {{1/6,-1/3,1/6},{(-33^0.5-1)/12,(-33^0.5+3)/18,(-33^0.5+15)/36},{(33^0.5-1)/12,(33^0.5+3)/18,(33^0.5+15)/36}}

  //var y = Expression.getFormaDeJordan(...);

</script>
```

Usage from node.js:
===================
```javascript
// npm install @yaffle/expression --save

  global.BigIntWrapper = require('@yaffle/expression/BigIntWrapper.js').BigIntWrapper;
  global.BigInteger = require('@yaffle/expression/BigInteger.js').BigInteger;
  require('@yaffle/expression/BigIntegerMath.js');
  global.Expression = require('@yaffle/expression/Expression.js').Expression;
  require('@yaffle/expression/toDecimalString.js');
  require('@yaffle/expression/GF2.js');
  require('@yaffle/expression/sin.js');
  require('@yaffle/expression/complex.js');
  global.NonSimplifiedExpression = require('@yaffle/expression/NonSimplifiedExpression.js').NonSimplifiedExpression;
  global.ExpressionParser = require('@yaffle/expression/ExpressionParser.js').ExpressionParser;
  global.Polynom = require('@yaffle/expression/Polynom.js').Polynom;
  global.Matrix = require('@yaffle/expression/Matrix.js').Matrix;
  require('@yaffle/expression/toLaTeX.js');
  require('@yaffle/expression/polynomial-roots-finding.js');
  require('@yaffle/expression/FormaCanonicalDeJordan.js');
  require('@yaffle/expression/ExpressionExtensions.js');
  require('@yaffle/expression/Condition.js');
  require('@yaffle/expression/toMathML.js');

  var p = Polynom.toPolynomial(ExpressionParser.parse("10x^5−17x^4−505x^3+1775x^2−249x−630"), ExpressionParser.parse("x"));
  console.log(p.getroots().toString()); // -1/2,5,21/5,(-73^0.5-7)/2,(73^0.5-7)/2

  var p = Polynom.toPolynomial(ExpressionParser.parse("x^5−2x^4−11x^3+26x^2−2x−13"), ExpressionParser.parse("x"));
  console.log(p.getZeros(5).result.toString()); // -3.412,-0.609,1.075,1.925,3.021

  var matrix = ExpressionParser.parse('{{1,2,3},{4,5,6},{7,8,9}}').matrix;
  console.log('matrix: ' + matrix.toString()); // matrix: {{1,2,3},{4,5,6},{7,8,9}}

  var x = Expression.getEigenvalues(matrix);
  var multiplicities = x.multiplicities;
  var eigenvalues = x.eigenvalues;
  console.log('eigenvalues: ' + eigenvalues.toString()); // eigenvalues: 0,(-3*33^0.5+15)/2,(3*33^0.5+15)/2

  var eigenvectors = Expression.getEigenvectors(matrix, x.eigenvalues).eigenvectors;
  console.log('eigenvectors: ' + eigenvectors.toString()); // eigenvectors: {{1},{-2},{1}},{{(-3*33^0.5-11)/22},{(-3*33^0.5+11)/44},{1}},{{(3*33^0.5-11)/22},{(3*33^0.5+11)/44},{1}}

  var y = Expression.diagonalize(matrix, eigenvalues, multiplicities, eigenvectors);
  console.log('diagonalization: ' + matrix.toString() + ' = ' + y.T.toString() + " * " + y.L.toString() + " * " + y.T_INVERSED.toString()); // diagonalization: {{1,2,3},{4,5,6},{7,8,9}} = {{1,(-3*33^0.5-11)/22,(3*33^0.5-11)/22},{-2,(-3*33^0.5+11)/44,(3*33^0.5+11)/44},{1,1,1}} * {{0,0,0},{0,(-3*33^0.5+15)/2,0},{0,0,(3*33^0.5+15)/2}} * {{1/6,-1/3,1/6},{(-33^0.5-1)/12,(-33^0.5+3)/18,(-33^0.5+15)/36},{(33^0.5-1)/12,(33^0.5+3)/18,(33^0.5+15)/36}}

  //var y = Expression.getFormaDeJordan(...);

```


Types
=====
```
  BigInteger
    BigInteger.BigInt(number)
    BigInteger.BigInt(string)
    BigInteger.toNumber(a)
    a.toString(radix)
    BigInteger.add(a, b)
    BigInteger.subtract(a, b)
    BigInteger.multiply(a, b)
    BigInteger.divide(a, b)
    BigInteger.remainder(a, b)
    BigInteger.exponentiate(a, n)
    BigInteger.unaryMinus(a)
    BigInteger.lessThan(a, b)
    BigInteger.nthRoot(a, n)
    BigInteger.gcd(a, b)
    BigInteger.primeFactor(a)
  Matrix
    rows()
    cols()
    e(row, column)
    isSquare()
    map(mapFunction)
    transpose()
    scale(x)
    multiply(b)
    add(b)
    subtract(b)
    augment(b)
    rowReduce(...)
    swapRows(...)
    toRowEchelon(...)
    determinant()
    rank()
    inverse()
    toString()
    pow(n)
    eql()
  Polynomial
    getDegree()
    getCoefficient()
  ExpressionParser
    parse(string, context)
  Expression
    .ZERO
    .ONE
    .TWO
    .TEN
    #add
    #subtract
    #multiply
    #divide
    #pow
    #equals
    #toString()
    #toMathML()
    #toLaTeX()
      Expression.Integer
        integer
      Expression.Symbol
        symbol
      Expression.NthRoot
        radicand
      Expression.Matrix
        matrix
      Expression.Polynomial
        polynomial
      Expression.Sin
        argument
      Expression.Cos
        argument
      Expression.Complex
        real
        imaginary
```
Similar projects
================
  https://coffeequate.readthedocs.io/en/latest/usage/
  http://algebrite.org - this page contains the list of "JavaScript Computer Algebra Systems"
  https://github.com/sloisel/numeric
