import Expression from './Expression.js';
import Polynomial from './Polynomial.js';

function deepFreeze(x) {
  if (Object.freeze != null) {
    if (typeof x !== "object") {
      Object.freeze(x);
      for (var key in x) {
        if (Object.prototype.hasOwnProperty.call(x, key)) {
          deepFreeze(x[key]);
        }
      }
    }
  }
}

function Condition(array) {
  deepFreeze(this);
  this.array = array;
}

Condition.NEZ = " != 0";
Condition.EQZ = " == 0";

Condition.prototype._and = function (operator, e) {
  if (operator !== Condition.NEZ && operator !== Condition.EQZ) {
    throw new TypeError();
  }
  if (e == undefined) {
    throw new RangeError();
  }
  if (this === Condition.FALSE) {
    return this;
  }

  var contains = function (array, operator, e) {
    for (var i = 0; i < array.length; i += 1) {
      if (array[i].operator === operator && array[i].expression.equals(e)) {
        return true;
      }
    }
    return false;
  };

  if (e instanceof Expression.GF2Value) {
    return this._and(operator, e.value === 0 ? Expression.ZERO : Expression.ONE);//?
  }

  //!new 2019-12-15:
  //!substitute:  x = 0, sin(x) != 0
  if (Expression.has(e, Expression.Sin) || Expression.has(e, Expression.Cos) || Expression.has(e, Expression.Exponentiation)) {
    var that = this;
    e = Expression._map(function (x) {
      if (x instanceof Expression.Function || x instanceof Expression.Exponentiation && (!(x.b instanceof Expression.Integer) || !(x.a instanceof Expression.Symbol))) {
        var arg = new Expression.Symbol('arg');
        var result = that.andZero(arg.subtract(x instanceof Expression.Exponentiation ? x.b : x.a));
        for (var i = 0; i < result.array.length; i += 1) {//TODO: fix
          var y = result.array[i];
          if (y.operator === Condition.EQZ) {
            var polynomial = Polynomial.toPolynomial(y.expression, arg);
            if (polynomial.getDegree() === 1) {
              var yy = polynomial.getCoefficient(0).negate().divide(polynomial.getCoefficient(1));
              if (!Expression.has(yy, Expression.Function)) {// sin(yy)/cos(yy) is supported
                if (x instanceof Expression.Exponentiation) {
                  //?TODO: tests
                  return x.a.pow(yy);
                } else {
                  yy = Expression.isConstant(yy) && !(yy.equals(Expression.ZERO)) ? new Expression.Radians(yy) : yy;
                  if (x instanceof Expression.Sin) {
                    return yy.sin();
                  } else if (x instanceof Expression.Cos) {
                    return yy.cos();
                  } else {
                    throw new TypeError("NotSupportedError");
                  }
                }
              }
            }
          }
        }
        return x;
      }
      return x;
    }, e);
  }
  //!

  var i = new Expression.Symbol('_i');
  function replaceSinCos(e) {
    return Expression._map(function (x) {
      if (x instanceof Expression.Sin) {
        var a = x.a;
        // Euler's formula
        return i.multiply(a).exp().subtract(i.multiply(a).negate().exp()).divide(Expression.TWO.multiply(Expression.I));
      }
      if (x instanceof Expression.Cos) {
        var a = x.a;
        // Euler's formula
        return i.multiply(a).exp().add(i.multiply(a).negate().exp()).divide(Expression.TWO);
      }
      return x;
    }, e);
  }
  function replaceBySinCos(e) {
    return Expression._map(function (x) {
      if (x instanceof Expression.Exponentiation && x.a === Expression.E) {
        var b = x.b;
        var p = Polynomial.toPolynomial(b, i);
        if (p.getDegree() === 1) {
          var q = p.getCoefficient(0);
          var w = p.getCoefficient(1);
          // Euler's formula
          return w.cos().add(Expression.I.multiply(w.sin())).multiply(q.exp());
        }
        if (p.getDegree() > 1) {
          throw new TypeError("!?");
        }
      }
      return x;
    }, e);
  }

  var add = function (oldArray, y) {
    //TODO: y is const
    if (contains(oldArray, y.operator, y.expression)) {//!TODO: it should work even without this (?)
      return oldArray;
    }
    if (contains(oldArray, y.operator === Condition.EQZ ? Condition.NEZ : Condition.EQZ, y.expression)) {
      return null;
    }

    var operator = null;// to not use a variable from scope accidently
    var e = y.expression;//!

    //!new
    if (e.isNegative()) {
      return add(oldArray, {expression: e.negate(), operator: y.operator});
    }

    // (x-1)^(1/2)
    if (e instanceof Expression.Exponentiation// &&
        //e.b.getNumerator() instanceof Expression.Integer &&
        //!e.b.getDenominator().equals(Expression.ONE)
        ) {
      return add(oldArray, {expression: e.a, operator: y.operator});
    }

    // (4*k+1)^(1/2)+1
    if (e instanceof Expression.Addition &&
        e.a instanceof Expression.Exponentiation &&
        Expression.isConstant(e.b) && //!
        e.a.b.getDenominator() instanceof Expression.Integer &&
        !e.a.b.getDenominator().equals(Expression.ONE)) {
      if (e.a.b.getDenominator().remainder(Expression.TWO).equals(Expression.ZERO) && !e.b.isNegative()) {
        return add(oldArray, {expression: Expression.ONE, operator: y.operator});
      }
      //if (!e.b.negate().pow(e.a.b.inverse()).pow(e.a.b).equals(e.b.negate())) {
      //  return add(oldArray, {expression: Expression.ONE, operator: y.operator});
      //}
      //TODO: fix
      return add(oldArray, {expression: e.a.a.pow(e.a.b.getNumerator()).subtract(e.b.negate().pow(e.a.b.getDenominator())), operator: y.operator});
    }

    if (y.expression instanceof Expression.Multiplication && y.expression.b instanceof Expression.IdentityMatrix) {
      return add(oldArray, {expression: y.expression.a, operator: y.operator});
    }
    if (y.expression instanceof Expression.Division) {
      var tmp = oldArray;
      tmp = add(tmp, {expression: y.expression.a, operator: y.operator});
      if (tmp == null) {
        return null;
      }
      tmp = add(tmp, {expression: y.expression.b, operator: Condition.NEZ});
      return tmp;
    }
    /*if (y.expression instanceof Expression.Division && Expression.isConstant(y.expression.b)) {
      y = {
        expression: y.expression.a,
        operator: y.operator
      };
    }*/
    if (y.expression instanceof Expression.Integer || y.expression instanceof Expression.Complex) {
      if (y.operator === Condition.NEZ && y.expression.equals(Expression.ZERO) ||
          y.operator === Condition.EQZ && !y.expression.equals(Expression.ZERO)) {
        return null;
      }
      return oldArray;
    }
    //TODO: check code coverage, remove extra branches
    if (Expression.isConstant(y.expression) && !y.expression.equals(Expression.ZERO)) {
      if (y.operator === Condition.NEZ) {
        return oldArray;
      }
      if (y.operator === Condition.EQZ) {
        return null;
      }
    }
    if (y.expression instanceof Expression.NthRoot) {
      return add(oldArray, {expression: y.expression.a, operator: y.operator});
    }
    if (y.expression instanceof Expression.Multiplication) {
      if (y.operator === Condition.EQZ) {
        if (y.expression.a instanceof Expression.Integer && !y.expression.a.equals(Expression.ZERO)) {
          //TODO: fix - ?
          y = {
            expression: y.expression.b,
            operator: y.operator
          };
          return add(oldArray, y);
        }
      }
      if (y.operator === Condition.NEZ) {
        var tmp = oldArray;
        tmp = add(tmp, {expression: y.expression.a, operator: Condition.NEZ});
        if (tmp == null) {
          return null;
        }
        tmp = add(tmp, {expression: y.expression.b, operator: Condition.NEZ});
        return tmp;
      }
    }

    var p = Expression.getMultivariatePolynomial(y.expression);
    if (p != null) {
      var content = p.p.getContent();
      if (!content.equals(Expression.ONE) && !content.equals(Expression.ONE.negate())) {
        // content * y.expression.divide(content)
        if (y.operator === Condition.NEZ) {
          var tmp = add(oldArray, {expression: y.expression.divide(content), operator: Condition.NEZ});
          if (tmp == null) {
            return null;
          }
          return add(tmp, {expression: content, operator: Condition.NEZ});
        }
        while (p != null) {
          if (y.operator === Condition.EQZ) {
            var sf = p.p.getSquareFreePolynomial();
            if (sf.getDegree() !== p.p.getDegree()) {
              //?
              return add(oldArray, {expression: y.expression.divide(p.p.divideAndRemainder(sf).quotient.calcAt(p.v)), operator: Condition.EQZ});
            }
          }
          content = p.p.getContent();
          p = Expression.getMultivariatePolynomial(content);
        }
        //!new
        if (Expression.isConstant(content)) {
          if (!content.equals(Expression.ONE) && !content.equals(Expression.ONE.negate())) {
            return add(oldArray, {expression: y.expression.divide(content), operator: Condition.EQZ});
          }
        }
        //y = {
        //  expression: y.expression.divide(content),
        //  operator: y.operator
        //};
      }
      if (p != null && p.p.getDegree() > 1 && p.p.getCoefficient(0).equals(Expression.ZERO)) {
        if (y.operator === Condition.NEZ) {
          var tmp = add(oldArray, {expression: p.v, operator: Condition.NEZ});
          if (tmp == null) {
            return null;
          }
          return add(tmp, {expression: y.expression.divide(p.v), operator: Condition.NEZ});
        }
      }
    }

    //!new 2019-12-24:
    var p = Expression.getMultivariatePolynomial(y.expression);
    if (p != null && p.p.getDegree() > 1) {
      var sf = p.p.getSquareFreePolynomial();
      if (sf.getDegree() !== p.p.getDegree()) {//TODO: test, fix
        return add(oldArray, {expression: sf.calcAt(p.v), operator: y.operator});
      }
    }
    //!

    var addRest = function (newArray, oldArray, i, other) {
      for (var j = i + 1; j < oldArray.length; j += 1) {
        newArray = add(newArray, oldArray[j]);
        if (newArray == null) {
          return null;
        }
      }
      newArray = add(newArray, other);
      return newArray;
    };


    var newArray = [];
    for (var i = 0; i < oldArray.length; i += 1) {
      var x = oldArray[i]; // TODO: const
      
      
      // (e**(tx)-e**(-tx))/(2i)
      // (e**(tx)+e**(-tx))/2

      // sin(x)=0, cos(x)=0
      //!new 2020-01-01:
      if (Expression.has(x.expression, Expression.Sin) || Expression.has(x.expression, Expression.Cos)) {
        if (Expression.has(y.expression, Expression.Sin) || Expression.has(y.expression, Expression.Cos)) {
          var xx = {
            operator: x.operator,
            expression: replaceSinCos(x.expression)
          };
          var yy = {
            operator: y.operator,
            expression: replaceSinCos(y.expression)
          };
          var tmp = add(add([], xx), yy);
          if (tmp == null) {
            return null;
          }
          if (tmp.length === 1) {
            return addRest(newArray, oldArray, i, {operator: tmp[0].operator, expression: replaceBySinCos(tmp[0])});
          }
          //TODO: ?
          //for (var i = 0; i < tmp.length; i++) {
          //  newArray = add(newArray, {operator: tmp[1].operator, expression: replaceBySinCos(tmp[i])});
          //}
          //return addRest(newArray, oldArray, i, {operator: Condition.EQZ, expression: Expression.ZERO});
        }
      }
      //!
      
      if ((x.operator === Condition.NEZ && y.operator === Condition.EQZ ||
           x.operator === Condition.EQZ && y.operator === Condition.NEZ) &&
           (Expression.isSingleVariablePolynomial(x.expression.multiply(y.expression)) || true)) {
        var g = x.expression.gcd(y.expression);
        while (!g.equals(Expression.ONE) && !g.equals(Expression.ONE.negate())) {
          if (x.operator === Condition.EQZ) {
            x = {
              operator: x.operator,
              expression: x.expression.divide(g)
            };
            // the change may affect all previous conditions:
          } else { // y.operator === Condition.EQZ
            y = {
              operator: y.operator,
              expression: y.expression.divide(g)
            };
            // we have not checked the y agains the branches in the beginning of the "add"
          }
          newArray = add(newArray, x);
          if (newArray == null) {
            return null;
          }
          return addRest(newArray, oldArray, i, y);
          //g = x.expression.gcd(y.expression);
        }
        //if (x.operator === Condition.EQZ) {
        //  var tmp = y;
        //  y = x;
        //  x = tmp;
        //}
        if (x.operator === Condition.EQZ && Expression.isConstant(x.expression)) {
          return null;
        }
        if (y.operator === Condition.EQZ && Expression.isConstant(y.expression)) {
          return null;
        }
        //if (!Expression.isSingleVariablePolynomial(x.expression.multiply(y.expression))) {
        //  if (x.operator === Condition.EQZ) {
        //    newArray.push(y);
        //  } else {
        //    newArray.push(x);
        //  }
        //}
      }
      if (x.operator === Condition.NEZ && y.operator === Condition.EQZ && Expression.isSingleVariablePolynomial(x.expression.multiply(y.expression))) {
        y = y;
      } else if (x.operator === Condition.EQZ && y.operator === Condition.NEZ && Expression.isSingleVariablePolynomial(x.expression.multiply(y.expression))) {
        y = x;
      } else if (x.operator === Condition.EQZ && y.operator === Condition.EQZ && Expression.isSingleVariablePolynomial(x.expression.multiply(y.expression))) {
        var g = x.expression.gcd(y.expression);
        if (g instanceof Expression.Integer) {
          return null;
        }
        y = {
          operator: y.operator,
          expression: g
        };
        return addRest(newArray, oldArray, i, y);
      } else if (x.operator === Condition.NEZ && y.operator === Condition.NEZ && Expression.isSingleVariablePolynomial(x.expression.multiply(y.expression))) {
        var g = x.expression.gcd(y.expression);
        x = {
          operator: x.operator,
          expression: x.expression.divide(g)
        };
        if (!Expression.isConstant(x.expression)) {
          newArray.push(x);
        }
      } else { // !isSingleVariablePolynomial
        // TODO: use Expression.isSingleVariablePolynomial(x.expression.multiply(y.expression))) here, and remove in the branches above

        var p = null;
        var pOperator = null;
        var pp = null;
        var other = null;
        var px = Expression.getMultivariatePolynomial(x.expression);
        var py = Expression.getMultivariatePolynomial(y.expression);
        //var xy = Expression.getMultivariatePolynomial(x.expression.multiply(y.expression));

        //console.assert(px != null && py != null);

        if (//xy != null &&
            //x.operator === Condition.EQZ &&
            //y.operator === Condition.EQZ &&
            px != null && py != null) {
              
          
          //!new 2019-24-12
          if (px != null && py != null && px.v.equals(py.v) && px.p.getDegree() !== 1 && py.p.getDegree() !== 1 && x.operator === Condition.EQZ && y.operator === Condition.EQZ) {
            //TODO: test, fix
            var tmp1 = py.p.getDegree() >= px.p.getDegree() ? Polynomial.pseudoRemainder(py.p, px.p) : py.p;
            var tmp2 = tmp1.calcAt(px.v);
            var tmp = {expression: tmp2, operator: y.operator};
            newArray = add(newArray, tmp);
            if (newArray == null) {
              return null;
            }
            return addRest(newArray, oldArray, i, x);
          }
          //!?

          //if (px != null && px.p.getDegree() !== 1 && py == null) {
            //py = {p: Polynomial.toPolynomial(y.expression, px.v), v: px.v};
            //if (xy.v === py.v) {
            //  py = null;
            //}
          //}
          //if (px == null && py != null && py.p.getDegree() !== 1) {
            //px = {p: Polynomial.toPolynomial(x.expression, py.v), v: py.v};
            //if (xy.v === px.v) {
            //  px = null;
            //}
          //}
          //if (px == null && py == null) {//?TODO:
          //  px = {p: Polynomial.toPolynomial(x.expression, xy.v), v: xy.v};
          //  py = {p: Polynomial.toPolynomial(y.expression, xy.v), v: xy.v};
          //}

          if (y.operator === Condition.EQZ && py != null && py.p.getDegree() === 1 &&
              x.operator === Condition.EQZ && px != null && px.p.getDegree() === 1) {
            //TODO: fix !!!
            //TODO: test linear systems
            if (Polynomial.toPolynomial(y.expression, px.v).getDegree() === 0) {
              px = null;
            }
            if (Polynomial.toPolynomial(x.expression, py.v).getDegree() === 0) {
              py = null;
            }

            if (px != null && py != null) { // TODO: ?
              if (!(px.p.getCoefficient(1) instanceof Expression.Integer)) {
                px = null;
              }
              if (!(py.p.getCoefficient(1) instanceof Expression.Integer)) {
                py = null;
              }
            }

            if (px != null && py != null) {
            
            //if (px.v.symbol < py.v.symbol) {//!
            if (px.v.compare4Addition(py.v) < 0) {
              px = null;
            }
            //}
            
            }
          }

          if (y.operator === Condition.EQZ && py != null && py.p.getDegree() === 1) {
            pp = py;
            p = x.expression;
            pOperator = x.operator;
            other = y;
          }
          if (x.operator === Condition.EQZ && px != null && px.p.getDegree() === 1) {
            pp = px;
            p = y.expression;
            pOperator = y.operator;
            other = x;
          }
        }
        if (pp != null) {
          var polynomial = Polynomial.toPolynomial(p, pp.v);
          var a = polynomial.calcAt(pp.p.getCoefficient(0).negate().divide(pp.p.getCoefficient(1)));
          if (!a.equals(p) &&
              (pp.p.getCoefficient(1) instanceof Expression.Integer || polynomial.divideAndRemainder(pp.p, "undefined") != undefined)) {
            var tmp = {
              operator: pOperator,
              expression: a
            };
            newArray = add(newArray, tmp);
            if (newArray == null) {
              return null;
            }
            if (true) {
              return addRest(newArray, oldArray, i, other);
            } else {
              y = other;
            }
          } else {
            newArray.push(x);
          }
        } else {
          newArray.push(x);
        }
      }
    }
    newArray.push(y);
    return newArray;
  };
  var newArray = add(this.array, {
    operator: operator,
    expression: e
  });
  if (newArray == null) {
    return Condition.FALSE;
  }
  if (newArray.length === 0) {
    return Condition.TRUE;
  }
  
  return new Condition(newArray);
};

Condition.prototype.andNotZero = function (e) {
  return this._and(Condition.NEZ, e);
};
Condition.prototype.andZero = function (e) {
  return this._and(Condition.EQZ, e);
};
Condition.prototype.isFalse = function () {
  return this === Condition.FALSE;
};
Condition.prototype.isTrue = function () {
  return this === Condition.TRUE;
};
Condition.prototype.toString = function (options) {
  if (this === Condition.TRUE || this === Condition.FALSE) {
    // 1) no need; 2) no way to distinguish TRUE and FALSE
    throw new TypeError();
  }
  if (this.array.length === 0) {
    // assertion
    throw new TypeError();
  }
  var s = '';
  for (var i = 0; i < this.array.length; i += 1) {
    s += (i !== 0 ? ', ' : '') + this.array[i].expression.toString(options) + (this.array[i].operator === Condition.NEZ ? ' != 0' : '') + (this.array[i].operator === Condition.EQZ ? ' == 0' : '');
  }
  return s;
};

Condition.TRUE = new Condition(new Array(0));
Condition.FALSE = new Condition(undefined);

export default Condition;
