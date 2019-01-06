/*global Expression */

(function () {
"use strict";

  var Integer = Expression.Integer;
  var Symbol = Expression.Symbol;
  var Addition = Expression.Addition;
  var Multiplication = Expression.Multiplication;
  var Division = Expression.Division;
  var Exponentiation = Expression.Exponentiation;
  var BinaryOperation = Expression.BinaryOperation;

var separateSinCos = function (e) {
  if (!(e instanceof Multiplication)) {
    throw new Error();
  }
  var sinCos = undefined;
  var other = undefined;
  for (var multiplications = e.factors(); multiplications != undefined; multiplications = multiplications.next()) {
    var v = multiplications.value();
    if (v instanceof Sin || v instanceof Cos || 
        (v instanceof Exponentiation && (v.a instanceof Sin || v.a instanceof Cos))) {
      sinCos = sinCos == undefined ? v : sinCos.multiply(v);
    } else {
      other = other == undefined ? v : other.multiply(v);
    }
  }
  return {
    sinCos: sinCos == undefined ? Expression.ONE : sinCos,
    other: other == undefined ? Expression.ONE : other
  };
};

var expandMainOp = function (u) {
  return u;
};

var contractTrigonometryInternal = function (a, b) {
  // sin(a) * sin(b) = (cos(a - b) - cos(a + b)) / 2
  // sin(a) * cos(b) = (sin(a + b) + sin(a - b)) / 2
  // cos(a) * sin(b) = (sin(a + b) - sin(a - b)) / 2
  // cos(a) * cos(b) = (cos(a - b) + cos(a + b)) / 2
  var ax = a.a;
  var bx = b.a;
  if (a instanceof Sin && b instanceof Sin) {
    return ax.subtract(bx).cos().divide(Expression.TWO).subtract(ax.add(bx).cos().divide(Expression.TWO));
  }
  if (a instanceof Sin && b instanceof Cos) {
    return ax.add(bx).sin().divide(Expression.TWO).add(ax.subtract(bx).sin().divide(Expression.TWO));
  }
  if (a instanceof Cos && b instanceof Sin) {
    return ax.add(bx).sin().divide(Expression.TWO).subtract(ax.subtract(bx).sin().divide(Expression.TWO));
  }
  if (a instanceof Cos && b instanceof Cos) {
    return ax.subtract(bx).cos().divide(Expression.TWO).add(ax.add(bx).cos().divide(Expression.TWO));
  }
  throw new Error();
};

// page 306
var contractTrigonometryPower = function (u) {
  var b = u.a;
  if (!(b instanceof Sin) && !(b instanceof Cos)) {
    return u;
  }
  var e = contractTrigonometryInternal(b, b).multiply(u.divide(b.multiply(b)));
  return contractTrigonometryRules(e.getNumerator()).divide(e.getDenominator());
};

// page 318
var contractTrigonometryProduct = function (u) {
  var i = u.factors();
  var a = i.value();
  i = i.next();
  var b = i == undefined ? undefined : i.value();
  i = i == undefined ? undefined : i.next();
  var rest = i == undefined || i.e == undefined ? Expression.ONE : i.e;//TODO: fix
  
  if (a instanceof Exponentiation) {
    a = contractTrigonometryPower(a);
    return contractTrigonometryRules(a.multiply(b).multiply(rest));
  }
  if (b instanceof Exponentiation) {
    b = contractTrigonometryPower(b);
    return contractTrigonometryRules(a.multiply(b).multiply(rest));
  }
  // (a instanceof Sin || a instanceof Cos) && (b instanceof Sin || b instanceof Cos)
  var c = contractTrigonometryInternal(a, b);

  return contractTrigonometryRules(c.multiply(rest));
};

// page 317
var contractTrigonometryRules = function (u) {
  var v = expandMainOp(u);
  if (v instanceof Exponentiation) {
    return contractTrigonometryPower(v);
  }
  if (v instanceof Multiplication) {
    var tmp = separateSinCos(v);
    var c = tmp.other;
    var d = tmp.sinCos;
    if (d.equals(Expression.ONE)) {
      return v;
    }
    if (d instanceof Sin || d instanceof Cos) {
      return v;
    }
    if (d instanceof Exponentiation) {
      return expandMainOp(c.multiply(contractTrigonometryPower(d)));
    }
    if (d instanceof Multiplication) {
      return expandMainOp(c.multiply(contractTrigonometryProduct(d)));
    }
    throw new Error();
  }
  if (v instanceof Addition) {
    var s = Expression.ZERO;
    for (var additions = v.summands(); additions != undefined; additions = additions.next()) {
      var y = additions.value();
      if (y instanceof Multiplication || y instanceof Exponentiation) {
        s = s.add(contractTrigonometryRules(y));
      } else {
        s = s.add(y);
      }
    }
    return s;
  }
  return v;
};

var map = function (f, u) {
  if (u instanceof Integer) {
    return f(u);
  }
  if (u instanceof Addition) {
    return f(u.a).add(f(u.b));
  }
  if (u instanceof Multiplication) {
    return f(u.a).multiply(f(u.b));
  }
  if (u instanceof Division) {
    return f(u.a).divide(f(u.b));
  }
  if (u instanceof Exponentiation) {
    return f(u.a).pow(f(u.b));
  }
  if (u instanceof Sin) {
    return f(u.a).sin();
  }
  if (u instanceof Cos) {
    return f(u.a).cos();
  }
  if (u instanceof Expression.Matrix) {
    return new Expression.Matrix(u.matrix.map(function (e, i, j) {
      return map(f, e);
    }));
  }
  if (u instanceof Expression.Polynomial) {//TODO: test case
    return new Expression.Polynomial(u.polynomial.map(function (c, d) {
      return map(f, c);
    }));
  }
  if (u instanceof Expression.GF2Value) {
    return u;
  }
  if (u instanceof Expression.NthRoot) {
    return u;
  }
  if (u instanceof Expression.Negation) {
    return u;//?
  }
  if (u instanceof Expression.Complex) {
    return u;//?
  }
  throw new Error();
};

// page 303

var expandTrigonometryRulesInternal = function (a, b, type) {
  if (type === "cos") {
    // cos(a + b) = cos(a) * cos(b) - sin(a) * sin(b)
    return expandTrigonometryRules(a, "cos").multiply(expandTrigonometryRules(b, "cos")).subtract(expandTrigonometryRules(a, "sin").multiply(expandTrigonometryRules(b, "sin")));
  }
  if (type === "sin") {
    // sin(a + b) = sin(a) * cos(b) + cos(a) * sin(b)
    return expandTrigonometryRules(a, "sin").multiply(expandTrigonometryRules(b, "cos")).add(expandTrigonometryRules(a, "cos").multiply(expandTrigonometryRules(b, "sin")));
  }
  throw new Error(type);
};

var expandTrigonometryRules = function (A, type) {
  if (A instanceof Addition) {
    return expandTrigonometryRulesInternal(A.a, A.b, type);
  } else if (A instanceof Multiplication) {
    var a = A.a;
    var b = A.b;
    if (!(a instanceof Integer)) {
      throw new Error();
    }
    if (a.compareTo(Expression.ONE.negate()) === 0) {
      if (type === "cos") {
        return expandTrigonometryRules(b, type);
      }
      if (type === "sin") {
        return expandTrigonometryRules(b, type).negate();
      }
    }
    var c = a.compareTo(Expression.ZERO) > 0 ? Expression.ONE : Expression.ONE.negate();
    return expandTrigonometryRulesInternal(c.multiply(b), a.subtract(c).multiply(b), type);
  } else if (A instanceof Division) {
    var t = simplifyConstantValue(A, type);
    if (t != null) {
      return t;
    }
    var a = A.a;
    var b = A.b;
    if (a instanceof Addition) {
      return expandTrigonometryRulesInternal(a.a.divide(b), a.b.divide(b), type);
    }
  }
  if (A instanceof Symbol) {
    if (type === "cos") {
      return A.cos();
    }
    if (type === "sin") {
      return A.sin();
    }
  }
  throw new Error();
};

// CA and SC, EA, p. 303

var expandTrigonometry = function (u) {
  if (u instanceof Integer || u instanceof Symbol) {
    return u;
  }
  var v = map(expandTrigonometry, u);
  if (v instanceof Sin) {
    return expandTrigonometryRules(v.a, "sin");
  }
  if (v instanceof Cos) {
    return expandTrigonometryRules(v.a, "cos");
  }
  return v;
};

var contractTrigonometry = function (u) {
  if (u instanceof Integer || u instanceof Symbol) {
    return u;
  }
  var v = map(contractTrigonometry, u);
  if (v instanceof Division) {//
    return contractTrigonometry(v.getNumerator()).divide(v.getDenominator());
  }
  if (v instanceof Multiplication || v instanceof Exponentiation || v instanceof Addition) {//! Addition - ?
    return contractTrigonometryRules(v);
  }
  if (v instanceof Cos || v instanceof Sin) {
    return v;
  }
  if (v instanceof Integer) {
    return v;
  }
  return v;//?
  //throw new Error();
};

// page 323

var hasTrigonometry = function (e) {//TODO: remove
  if (e instanceof BinaryOperation) {
    return hasTrigonometry(e.a) || hasTrigonometry(e.b);
  }
  return e instanceof Cos || e instanceof Sin;
};

var simplifyTrigonometry = function (u) {
  if (!hasTrigonometry(u)) {
    return u;
  }
  var n = expandTrigonometry(u.getNumerator());
  n = contractTrigonometry(n);
  var d = expandTrigonometry(u.getDenominator());
  d = contractTrigonometry(d);
  return n.divide(d);
};

Expression.simplifyTrigonometry = simplifyTrigonometry;//?

Expression.Multiplication.prototype.compare4Multiplication = function (y) {
  var x = this;
  if (y instanceof Addition) {//TODO: fix
    return 0 - y.compare4Multiplication(x);
  }
  var i = x.factors();
  var j = y.factors();
  while (i != undefined && j != undefined) {
    var a = i.value();
    i = i.next();
    var b = j.value();
    j = j.next();
    var c = a.compare4Multiplication(b);
    if (c !== 0) {
      return c;
    }
  }
  return i != undefined ? +1 : (j != undefined ? -1 : 0);
};

Expression.Multiplication.compare4Addition = function (x, y) {
  var i = x.factors();
  var j = y.factors();
  while (i != undefined && j != undefined) {
    var a = i.value();
    i = i.next();
    var b = j.value();
    j = j.next();
    var c = a.compare4Addition(b);
    if (c !== 0) {
      return c;
    }
  }
  return i != undefined ? +1 : (j != undefined ? -1 : 0);
};


// cos(2 * x) * cos(x)
Expression.Multiplication.prototype.compare4MultiplicationSymbol = function (x) {
  return 0 - this.compare4Multiplication(x);
};

Expression.Function.prototype.compare4Addition = function (y) {
  if (y instanceof Expression.Function) {
    return this.name < y.name ? -1 : (y.name < this.name ? +1 : this.a.compare4Addition(y.a));
  }
  return +1;
};

Expression.prototype.compare4AdditionSymbol = function (x) {
  var y = this;
  if (y instanceof Expression.Function) {
    return -1;
  }
  return Expression.prototype.compare4Addition.call(x, y);
};

Symbol.prototype.compare4Addition = function (y) {
  return y.compare4AdditionSymbol(this);
};

Expression.Function.prototype.compare4Multiplication = function (y) {
  if (y instanceof Expression.NthRoot) {
    return +1;
  }
  if (y instanceof Expression.Function) {
    return this.name < y.name ? -1 : (y.name < this.name ? +1 : this.a.compare4Multiplication(y.a));
  }
  return +1;
};

Expression.Function.prototype.compare4MultiplicationInteger = function (x) {
  return 0 - this.compare4Multiplication(x);
};

Expression.Function.prototype.compare4MultiplicationSymbol = function (x) {
  return -1;//?
};

Expression.Function.prototype.compare4MultiplicationNthRoot = function (x) {
  return 0 - this.compare4Multiplication(x);
};

Expression.Function.prototype.pow = function (y) {
  if (this instanceof Expression.NthRoot) {
    return Expression.prototype.pow.call(this, y);
  }
  if (y instanceof Expression.Integer) {
    if (y.compareTo(Expression.ONE) > 0) {
      return new Expression.Exponentiation(this, y);
    }
    return Expression.prototype.pow.call(this, y);
  }
  throw new RangeError("NotSupportedError");
};

function Sin(x) {
  Expression.Function.call(this, "sin", x);
}
Sin.prototype = Object.create(Expression.Function.prototype);

//TODO: new 2017-04-26
var simplifyConstantValueInternal = function (a, type) {
  a = a.remainder(Integer.parseInteger("360"));
  var d = Number.parseInt(a.toString(), 10);
  if (type === "cos") {
    d = 90 - d;
    if (d >= 360 - 90) {
      d -= 360;
    }
  }
  if (d < 0) {
    d += 360;
  }
  var s = d >= 180 ? Expression.ONE.negate() : Expression.ONE;
  if (d >= 180) {
    d -= 180;
  }
  if (d >= 90) {
    d = 180 - d;
  }
  if (d === 0) {
    return s.multiply(Expression.ZERO);
  }
  if (d === 15) {
    return s.multiply(Expression.TWO.add(Expression.ONE).squareRoot().subtract(Expression.ONE).divide(Expression.TWO.multiply(Expression.TWO.squareRoot())));
  }
  if (d === 30) {
    return s.multiply(Expression.ONE.divide(Expression.TWO));
  }
  if (d === 45) {
    return s.multiply(Expression.ONE.divide(Expression.TWO.squareRoot()));
  }
  if (d === 60) {
    return s.multiply(Expression.ONE.add(Expression.TWO).squareRoot().divide(Expression.TWO));
  }
  if (d === 75) {
    return s.multiply(Expression.TWO.add(Expression.ONE).squareRoot().add(Expression.ONE).divide(Expression.TWO.multiply(Expression.TWO.squareRoot())));
  }
  if (d === 90) {
    return s.multiply(Expression.ONE);
  }
  if (d === 18) {
    return Expression.TWO.add(Expression.TWO).add(Expression.ONE).squareRoot().subtract(Expression.ONE).divide(Expression.TWO.add(Expression.TWO));
  }
  if (d === 54) {
    return Expression.TWO.add(Expression.TWO).add(Expression.ONE).squareRoot().add(Expression.ONE).divide(Expression.TWO.add(Expression.TWO));
  }
  return undefined;
};

var simplifyConstantValue = function (x, type) {
  if (x instanceof Integer) {
    if (x.compareTo(Expression.ZERO) === 0) {
      return simplifyConstantValueInternal(Expression.ZERO, type);
    }
  }
  if (x instanceof Expression.Degrees) {
    return simplifyConstantValueInternal(x.value.simplify(), type);
  }
  var a = undefined;
  var b = undefined;
  if (x === Expression.PI) {
    a = Expression.ONE;
    b = Expression.ONE;
  } else if (x instanceof Multiplication && x.a instanceof Integer && x.b === Expression.PI) {
    a = x.a;
    b = Expression.ONE;
  } else if (x instanceof Division && x.b instanceof Integer && x.a === Expression.PI) {
    a = Expression.ONE;
    b = x.b;
  } else if (x instanceof Division && x.b instanceof Integer && x.a instanceof Multiplication && x.a.a instanceof Integer && x.a.b === Expression.PI) {
    a = x.a.a;
    b = x.b;
  }
  if (a != undefined && b != undefined) {
    b = Number.parseInt(b.toString(), 10);
    if (b === 1 || b === 2 || b === 3 || b === 4 || b === 6 || b === 12) {
      var d = a.multiply(Integer.parseInteger(Math.floor(180 / b).toString()));
      return simplifyConstantValueInternal(d, type);
    }
  }
  return undefined;
};

var isArgumentValid = function (x) {
  if (simplifyConstantValue(x, "sin") != undefined) {
    return true;
  }
  if (!(Expression.isScalar(x) && x instanceof Symbol) &&
      !(x instanceof Multiplication && x.a instanceof Integer && Expression.isScalar(x.b) && x.b instanceof Symbol) &&
      !(x instanceof Addition && isArgumentValid(x.a) && isArgumentValid(x.b)) &&
      !(x instanceof Division && x.b instanceof Integer && x.a instanceof Addition && isArgumentValid(x.a.a.divide(x.b)) && isArgumentValid(x.a.b.divide(x.b)))) {
    return false;
  }
  return true;
};

Expression.prototype.sin = function () {
  var x = this;
  if (!isArgumentValid(x)) {
    throw new RangeError("NotSupportedError");
  }
  if (x.isNegative()) {
    return new Sin(x.negate()).negate();
  }
  var t = simplifyConstantValue(x, "sin");
  if (t != undefined) {
    return t;
  }
  return new Sin(x);
};

function Cos(x) {
  Expression.Function.call(this, "cos", x);
}
Cos.prototype = Object.create(Expression.Function.prototype);

Expression.prototype.cos = function () {
  var x = this;
  if (!isArgumentValid(x)) {
    throw new RangeError("NotSupportedError");
  }
  if (x.isNegative()) {
    return new Cos(x.negate());
  }
  var t = simplifyConstantValue(x, "cos");
  if (t != undefined) {
    return t;
  }
  return new Cos(x);
};

Expression.simplifications.push(simplifyTrigonometry);

Expression.Sin = Sin;
Expression.Cos = Cos;

//Expression.Negation.prototype.compare4Multiplication = function (y) {
//TODO: fix, more tests
//  return new Expression.Multiplication(Expression.ONE.negate(), this.a).compare4Multiplication(y);
//};

Expression.Addition.prototype.compare4Addition = function (y) {
  // cos(a + b) + cos(a + b)
  var x = this;
  return Expression.Addition.compare4Addition(x, y);
};

//!!!
Expression.Addition.prototype.compare4Multiplication = function (y) {
  var x = this;
  var i = x.summands();
  var j = y.summands();
  while (i != undefined && j != undefined) {
    var a = i.value();
    i = i.next();
    var b = j.value();
    j = j.next();
    var c = a.compare4Multiplication(b);
    if (c !== 0) {
      return c;
    }
  }
  return i != undefined ? +1 : (j != undefined ? -1 : 0);
};

Expression.Addition.prototype.compare4MultiplicationSymbol = function (x) {
  return 0 - this.compare4Multiplication(x);
};

Expression.Addition.compare4Addition = function (x, y) {
  var i = x.summands();
  var j = y.summands();
  while (i != undefined && j != undefined) {
    var a = i.value();
    i = i.next();
    var b = j.value();
    j = j.next();
    var c = a.compare4Addition(b);
    if (c !== 0) {
      return c;
    }
  }
  return i != undefined ? +1 : (j != undefined ? -1 : 0);
};

//!!!


//!new 2017-04-26
Expression.Degrees = function (value) {
  this.value = value;
};
Expression.Degrees.prototype = Object.create(Expression.prototype);
Expression.Degrees.prototype.toString = function (options) {
  return this.value.toString(options) + "\u00B0";
};
Expression.Degrees.prototype.toMathML = function (options) {
  return "<mrow>" + this.value.toMathML(options) + "<mo>&#x2062;</mo><mi>&deg;</mi></mrow>";
};

}());
