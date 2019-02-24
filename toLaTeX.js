import Expression from './Expression.js';
import NonSimplifiedExpression from './Expression.js';

Expression.prototype.toLaTeX = function (options) {
  throw new RangeError();
};

Expression.Integer.prototype.toLaTeX = function (options) {
  return this.value.toString();
};
Expression.Function.prototype.toLaTeX = function (options) {
  return "\\" + this.name + "\\left(" + this.a.toLaTeX(options) + "\\right)";
};
Expression.SquareRoot.prototype.toLaTeX = function (options) {
  return "\\sqrt{" + this.a.toLaTeX(options) + "}";
};
Expression.NthRoot.prototype.toLaTeX = function (options) {
  return "\\sqrt[" + this.n + "]{" + this.a.toLaTeX(options) + "}";
};
Expression.BinaryOperation.prototype.toLaTeX = function (options) {
  var a = this.a;
  var b = this.b;
  var fa = a.getPrecedence() + (a.isRightToLeftAssociative() ? -1 : 0) < this.getPrecedence();
  var fb = this.getPrecedence() + (this.isRightToLeftAssociative() ? -1 : 0) >= b.getPrecedence();
  fa = fa || a.isUnaryPlusMinus();
  fb = fb || b.isUnaryPlusMinus(); // 1*-3 -> 1*(-3)
  fa = fa || (this instanceof Expression.Exponentiation && a instanceof Expression.Function); // cos(x)^(2+3)
  return (fa ? "\\left(" : "") + a.toLaTeX(options) + (fa ? "\\right)" : "") +
         this.getS() +
         (fb ? "\\left(" : "") + b.toLaTeX(options) + (fb ? "\\right)" : "");
};
Expression.Negation.prototype.toLaTeX = function (options) {
  var b = this.b;
  var fb = this.getPrecedence() + (this.isRightToLeftAssociative() ? -1 : 0) >= b.getPrecedence();
  fb = fb || b.isUnaryPlusMinus();
  // assert(fa === false);
  return "-" + (fb ? "\\left(" : "") + b.toLaTeX(options) + (fb ? "\\right)" : "");
};
Expression.Division.prototype.toLaTeX = function (options) {
  return "\\frac{" + this.a.toLaTeX(options) + "}{" + this.b.toLaTeX(options) + "}";
};
Expression.Symbol.prototype.toLaTeX = function (options) {
  return this.symbol;
};
Expression.Matrix.prototype.toLaTeX = function (options) {
  var x = this.matrix;
  var s = "";
  s += "\\begin{pmatrix}\n";
  for (var i = 0; i < x.rows(); i += 1) {
    for (var j = 0; j < x.cols(); j += 1) {
      var e = x.e(i, j);
      s += e.toLaTeX(options) + (j + 1 < x.cols() ? " & " : (i + 1 < x.rows() ? " \\\\" : "") + "\n");
    }
  }
  s += "\\end{pmatrix}";
  return s;
};
Expression.Complex.prototype.toLaTeX = function (options) {
  return this.toString(options);
};
NonSimplifiedExpression.prototype.toLaTeX = function (options) {
  return this.e.toLaTeX(options);
};
