
import toDecimalStringInternal from './toDecimalString.js';
import Polynomial from './Polynomial.js';//TODO:
import Expression from './Expression.js';//TODO:
import NonSimplifiedExpression from './NonSimplifiedExpression.js';//TODO:

//TODO: remove - ?
import Condition from './Condition.js';


//!TODO: remove
Polynomial.prototype.toString = function (options) {
  options = options || {};
  return this.toExpression(options.polynomialVariable || (new Expression.Symbol("x"))).toString(options);
};
Polynomial.prototype.toMathML = function (options) {
  options = options || {};
  return this.toExpression(options.polynomialVariable || (new Expression.Symbol("x"))).toMathML(options);
};
Expression.Polynomial.prototype.toString = function (options) {
  return this.polynomial.toString(options);
};
Expression.Polynomial.prototype.toMathML = function (options) {
  return this.polynomial.toMathML(options);
};

// coefficient - Expression
// variable - Expression
var printPartOfAddition = function (isLast, isFirst, coefficient, variable, options) {
  if (coefficient.equals(Expression.ZERO)) {
    return (isLast && isFirst ? Expression.ZERO.toMathML() : "");
  }
  var isNegative = false;
  if (coefficient.isNegative() && !isFirst) {
    isNegative = true;
    coefficient = coefficient.negateCarefully();//?
  }
  var precedenceOfMultiptication = new Expression.Multiplication(Expression.ZERO, Expression.ZERO).getPrecedence();
  var areBracketsRequired = coefficient.getPrecedence() < precedenceOfMultiptication; //?
  var c = coefficient.equals(Expression.ONE);
  //TODO: fix
  return (isFirst ? '' : '') +
         (!isFirst && isNegative ? '<mo form="infix">&minus;</mo>' : '') +
         (!isFirst && !isNegative ? '<mo form="infix">+</mo>' : '') +
         (c ? '' : '<mrow>') +
         (c || !areBracketsRequired ? '' : '<mrow><mo>(</mo>') +
         (c ? '' : coefficient.toMathML(options)) +
         (c || !areBracketsRequired ? '' : '<mo>)</mo></mrow>') +
         (c ? '' : '<mo>&sdot;</mo>') +
         variable.toMathML(options) +
         (c ? '' : '</mrow>');
};


var decimalToMathML = function (parts) {
  var formatNumber = function (s) {
    var locale = undefined;//TODO: ?
    var numberFormat = getNumberFormat(locale);
    if (numberFormat == undefined) {
      return s;
    }
    return numberFormat.format(s + '1').slice(0, -numberFormat.format('1').length);
  };
  return (parts.exponentInteger !== "" ? "<mrow>" : "") +
         (parts.plusSign !== "" || parts.minusSign !== "" ? "<mrow>" : "") +
         (parts.plusSign !== "" ? "<mo>+</mo>" : "") +
         (parts.minusSign !== "" ? "<mo>&minus;</mo>" : "") +
         "<mn>" + formatNumber(parts.number) + "</mn>" +
         (parts.plusSign !== "" || parts.minusSign !== "" ? "</mrow>" : "") +
         (parts.exponentInteger !== "" ? "<mo lspace=\"0\" rspace=\"0\">&sdot;</mo><msup>" + "<mn>" + formatNumber('10') + "</mn>" + decimalToMathML({plusSign: "", minusSign: parts.exponentMinusSign, number: parts.exponentInteger, exponentMinusSign: "", exponentInteger: ""}) + "</msup>" : "") +
         (parts.exponentInteger !== "" ? "</mrow>" : "");
};

var complexToMathML = function (real, imaginary) {
  if (imaginary.replace(/<[^>]+>/g, '') === '1') {
    return '<mrow>' + real + '<mo>+</mo><mi>&ii;</mi></mrow>';
  }
  if (imaginary.replace(/<[^>]+>/g, '').replace(/&minus;/g, '-') === '-1') {
    return '<mrow>' + real + '<mo>&minus;</mo><mi>&ii;</mi></mrow>';
  }
  if (real === '') {
    return '<mrow>' + imaginary + '<mo>' + (/<msup>/.test(imaginary) || /[^0-9]/.test(imaginary) ? '&sdot;' : '&it;') + '</mo><mi>&ii;</mi></mrow>';
  }
  var signBetween = '+';
  if (/<mrow><mo>&minus;<\/mo><mn>([^<]*)<\/mn><\/mrow>(?!<\/msup>)/.test(imaginary)) {
    signBetween = '-';
    imaginary = imaginary.replace(/<mrow><mo>&minus;<\/mo><mn>([^<]*)<\/mn><\/mrow>(?!<\/msup>)/g, '<mn>$1</mn>');
  }
  return '<mrow>' + real + '<mo>' + (signBetween === '-' ? '&minus;' : '+') + '</mo>' + '<mrow>' + imaginary + '<mo>' + (/<msup>/.test(imaginary) || /[^0-9]/.test(imaginary) ? '&sdot;' : '&it;') + '</mo>' + '<mi>&ii;</mi>' + '</mrow>' + '</mrow>';
};
Expression._decimalToMathML = decimalToMathML;
Expression._complexToMathML = complexToMathML;

  function isConstant(e) {
    if (e instanceof Expression.Symbol) {
      //return false;
      return e === Expression.E || e === Expression.PI;
    }
    if (e instanceof Expression.NonSimplifiedExpression) {
      return false;
    }
    if (e instanceof Expression.Matrix) {
      return false;
    }
    if (e instanceof Expression.Polynomial) {
      return false;
    }
    if (e instanceof Expression.BinaryOperation) {
      //return ((e.a === Expression.E || e.a === Expression.PI) || isConstant(e.a)) && isConstant(e.b);
      return isConstant(e.a) && isConstant(e.b);
    }
    if (e instanceof Expression.Negation) {
      return isConstant(e.b);
    }
    if (e instanceof Expression.Function) {
      return isConstant(e.a);
    }
    if (e instanceof Expression.Radians) {
      return isConstant(e.value);
    }
    return true;
  }

var groupByTerm = function (e) {
  if (e instanceof Expression.Division || e instanceof Expression.Addition) {
    var numerator = e.getNumerator();
    var denominator = e.getDenominator();
    if (denominator instanceof Expression.Integer) {
      var summands = [];
      for (var summand of numerator.summands()) {
        summands.push(summand);
      }
      var map = {};//TODO: FIX, test(?)
      summands.reverse();
      for (var i = 0; i < summands.length; i += 1) {
        var summand = summands[i];
        //var constant = Expression.getConstant(summand); - ?
        var constant = Expression.ONE;
        for (var factor of summand.factors()) {
          if (isConstant(factor)) {
            constant = constant.multiply(factor);
          }
        }
        var term = summand.divide(constant);
        var key = "_" + term.toString();
        map[key] = map[key] || {
          constant: Expression.ZERO,
          term: term
        };
        map[key].constant = map[key].constant.add(constant.divide(denominator));
      }
      var result = null;
      for (var key in map) {
        if (Object.prototype.hasOwnProperty.call(map, key)) {
          var constant = map[key].constant;
          var term = map[key].term;
          var s = false;
          if (result != null && constant.isNegative()) {
            constant = constant.negate();
            s = true;
          }
          var x = new NonSimplifiedExpression(term.equals(Expression.ONE) ? constant : (constant.equals(Expression.ONE) ? map[key].term : new Expression.Multiplication(constant, map[key].term)));
          result = result == null ? x : new Expression.Addition(result, s ? new Expression.Multiplication(Expression.ONE.negate(), x) : x);
        }
      }
      return new NonSimplifiedExpression(result);
    }
  }
  return undefined;
};

var getRounding = function (options) {
  var rounding = options != undefined ? options.rounding : undefined;
  if (options != undefined && rounding == undefined && options.fractionDigits != undefined && options.fractionDigits !== -1) {
    console.debug('deprecated fractionDigits, use rounding instead');
    rounding = {fractionDigits: options.fractionDigits};
  }
  return rounding;
};

//TODO: move
Expression.toDecimalString = function (x, options) {
  var rounding = getRounding(options);
  if (rounding != null) {
    if (isConstant(x)) {
      return toDecimalStringInternal(x, rounding, decimalToMathML, complexToMathML);
    } else if (!Expression.has(x, NonSimplifiedExpression)) {
      var grouped = groupByTerm(x);
      if (grouped != undefined) {
        return grouped.toMathML(options);
      }
    }
  }
  return undefined;
};

var getPrecedence = function (x, options) {
  var rounding = getRounding(options);
  if (rounding != null && isConstant(x.unwrap()) && Expression.has(x.unwrap(), Expression.Complex)) {
    if (!x.unwrap().equals(Expression.I)) {
      return new Expression.Addition(Expression.ONE, Expression.ONE).getPrecedence();
    }
  }
  return x.getPrecedence();
};

Expression.idCounter = 0;
Expression.id = function () {
  return (Expression.idCounter += 1).toString();
};

//TODO: ?
Expression.escapeHTML = function (s) {
  return s.replace(/&/g, "&amp;")
          .replace(/"/g, "&quot;")
          .replace(/</g, "&lt;")
          .replace(/>/g, "&gt;");
};

Expression.Matrix.prototype.toMathML = function (options) {
  var x = this.matrix;
  options = Expression.setTopLevel(true, options);

  var useMatrixContainer = options.useMatrixContainer == undefined ? true : options.useMatrixContainer;
  //TODO: fix!
  var braces = options.useBraces == undefined ? undefined : options.useBraces;
  var columnlines = options.columnlines == undefined ? 0 : options.columnlines;
  var variableNames = options.variableNames == undefined ? undefined : options.variableNames;

  var verticalStrike = options.verticalStrike == undefined ? -1 : options.verticalStrike;
  var horizontalStrike = options.horizontalStrike == undefined ? -1 : options.horizontalStrike;

  var cellIdGenerator = options.cellIdGenerator == undefined ? undefined : options.cellIdGenerator;
  var pivotCell = options.pivotCell == undefined ? undefined : options.pivotCell;

  var isLUDecomposition2 = options.isLUDecomposition2 == undefined ? undefined : options.isLUDecomposition2;
  var highlightRow = options.highlightRow == undefined ? -1 : options.highlightRow;
  var highlightCol = options.highlightCol == undefined ? -1 : options.highlightCol;

  options = Object.assign({}, options, {
    useBraces: undefined,
    columnlines: undefined,
    variableNames: undefined,
    verticalStrike: undefined,
    horizontalStrike: undefined,
    cellIdGenerator: undefined,
    pivotCell: undefined,
    isLUDecomposition2: undefined,
    highlightRow: undefined,
    highlightCol: undefined
  });

  var result = "";
  var rows = x.rows();
  var cols = x.cols();
  var i = -1;

  //TODO: remove `Expression.id()`
  var containerId = options.idPrefix + "-" + Expression.id();
  if (useMatrixContainer) {
    result += "<munder accentunder=\"true\">";
    // <menclose href="#"> will not be supported by MathML Core, so using <mrow>
    result += '<mrow id="' + containerId + '" data-matrix="' + Expression.escapeHTML(x.toString()) + '" draggable="true" tabindex="0" contextmenu="matrix-menu">';
  }

  result += braces == undefined ? '<mrow><mo>(</mo>' : '<mrow>' + (braces[0] === ' ' ? '' : '<mo>' + braces[0] + '</mo>');
  var columnlinesAttribute = "";
  if (columnlines !== 0 && cols - 1 > 0) {
    var k = -1;
    while (++k < cols - 1) {
      columnlinesAttribute += (cols - 1 + columnlines === k ? "solid " : "none ");
    }
    // whitespace
    columnlinesAttribute = columnlinesAttribute.slice(0, -1);
  }
  //! 2017-07-06 rowspacing="0ex" was added to make it look better with Native MathML (when it is supported) and to have the same style as in mathml.css
  //  rowspacing="0ex" is also needed when verticalStrike !== -1
  var useColumnspacing = verticalStrike !== -1 || horizontalStrike !== -1 || pivotCell != undefined || cellIdGenerator != undefined;
  result += "<mtable" + " rowspacing=\"0ex\"" + (useColumnspacing ? " columnspacing=\"0em\"" : "") + (variableNames != undefined ? " columnalign=\"right\"" : "") + (columnlinesAttribute !== "" ? " columnlines=\"" + columnlinesAttribute + "\"" : "") + ">";
  while (++i < rows) {
    var j = -1;
    if (variableNames != undefined) {// TODO: fix?
      //TODO: use code from polynomialToExpression3 (shared)
      var row = "";
      var wasNotZero = false;
      while (++j < cols - 1) {
        // TODO: fix `new Expression.Symbol()`
        row += "<mtd>";
        row += printPartOfAddition(j === cols - 2, !wasNotZero, x.e(i, j), new Expression.Symbol(variableNames[j]), options);
        row += "</mtd>";
        wasNotZero = wasNotZero || !x.e(i, j).equals(Expression.ZERO);
      }
      row += "<mtd><mo>=</mo></mtd><mtd>" + x.e(i, cols - 1).toMathML(options) + "</mtd>";
      if (wasNotZero || !x.e(i, cols - 1).equals(Expression.ZERO)) {
        result += "<mtr>";
        result += row;
        result += "</mtr>";
      }
    } else {
      result += "<mtr>";
      while (++j < cols) {
        result += "<mtd" + (cellIdGenerator != undefined ? " id=\"" + cellIdGenerator(i, j) + "\"" : "") + ">";
        if (pivotCell != undefined && i === pivotCell.i && j === pivotCell.j) {
          result += "<mstyle mathvariant=\"bold\">";
          result += "<menclose notation=\"circle\">";
        }
        if (horizontalStrike === i) {
          result += "<menclose notation=\"horizontalstrike\">";
        }
        if (verticalStrike === j) {
          result += "<menclose notation=\"verticalstrike\">";
        }
        if (useColumnspacing) {
          result += "<mpadded width=\"+0.8em\" lspace=\"+0.4em\">";
        }
        var highlight = j < i && isLUDecomposition2 ||
                        highlightRow === i && (columnlines === 0 || j <= cols - 1 + columnlines) || highlightCol === j;
        if (highlight) {
          result += "<mrow mathbackground=\"#80FF80\" mathcolor=\"#3C78C2\">";
        }
        result += x.e(i, j).toMathML(options);
        if (highlight) {
          result += "</mrow>";
        }
        if (useColumnspacing) {
          result += "</mpadded>";
        }
        if (verticalStrike === j) {
          result += "</menclose>";
        }
        if (horizontalStrike === i) {
          result += "</menclose>";
        }
        if (pivotCell != undefined && i === pivotCell.i && j === pivotCell.j) {
          result += "</menclose>";
          result += "</mstyle>";
        }
        result += "</mtd>";
      }
      result += "</mtr>";
    }
  }
  result += "</mtable>";
  result += braces == undefined ? '<mo>)</mo></mrow>' : (braces[1] === ' ' ? '' : '<mo>' + braces[1] + '</mo>') + '</mrow>';

  if (useMatrixContainer) {
    result += '</mrow>';

    result += "<mtext>";
    result += "<button type=\"button\" class=\"matrix-menu-show\" data-for-matrix=\"" + containerId + "\" aria-haspopup=\"true\"></button>";
    result += "</mtext>";

    result += "</munder>";
  }

  return result;
};

Expression.Determinant.prototype.toMathML = function (options) {
  var x = this;
  if (x.a instanceof Expression.Matrix || (x.a instanceof NonSimplifiedExpression && x.a.e instanceof Expression.Matrix)) {
    options = Object.assign({}, options, {
      useBraces: ["|", "|"]
    });
    //TODO: fix
    return x.a.toMathML(options);
  }
  return "<mrow><mo>|</mo>" + x.a.toMathML(options) + "<mo>|</mo></mrow>";
};
Expression.Transpose.prototype.toMathML = function (options) {
  var x = this;
  //TODO: ^T ?
  // https://www.w3.org/TR/MathML3/chapter4.html#contm.transpose
  var f = this.getPrecedence() >= x.a.getPrecedence();
  return "<msup>" +
         (f ? "<mrow><mo>(</mo>" : "") + x.a.toMathML(options) + (f ? "<mo>)</mo></mrow>" : "") +
         "<mi>T</mi>" +
         "</msup>";
};
Expression.SquareRoot.prototype.toMathML = function (options) {
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  return "<msqrt>" +
         this.a.toMathML(Expression.setTopLevel(true, options)) +
         "</msqrt>";
};
Expression.CubeRoot.prototype.toMathML = function (options) {
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  return "<mroot>" +
         this.a.toMathML(Expression.setTopLevel(true, options)) +
         "<mi>" + 3 + "</mi>" +
         "</mroot>";
};
Expression.NthRoot.prototype.toMathML = function (options) {
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  return "<mroot>" +
         this.a.toMathML(Expression.setTopLevel(true, options)) +
         "<mi>" + this.n + "</mi>" +
         "</mroot>";
};
Expression.Function.prototype.toMathML = function (options) {
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  var x = this;
  var fa = !(x.a instanceof Expression.Matrix) && !(x.a instanceof NonSimplifiedExpression && x.a.e instanceof Expression.Matrix);//?
  //TODO: fix
  return "<mrow>" +
         "<mi>" + (Expression.denotations[x.name] || x.name) + "</mi>" +
         "<mo>&af;</mo>" +
         (fa ? "<mrow><mo>(</mo>" : "") +
         x.a.toMathML(Expression.setTopLevel(true, options)) +
         (fa ? "<mo>)</mo></mrow>" : "") +
         "</mrow>";
};
Expression.Division.prototype.toMathML = function (options) {
  if (options != null && options.nofractions) {
    return Expression.BinaryOperation.prototype.toMathML.call(this, options);
  }
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  var x = this;
  var denominator = x.getDenominator();
  var numerator = x.getNumerator();
  //???
  //if (numerator.isNegative()) {
  //  return "<mrow><mo>&minus;</mo>" + x.negateCarefully().toMathML(options) + "</mrow>";
  //}
  return "<mfrac>" +
         numerator.toMathML(Expression.setTopLevel(true, options)) +
         denominator.toMathML(Expression.setTopLevel(true, options)) +
         "</mfrac>";
};

Expression.locale = undefined;
var cacheLocale = '';
var cachedNumberFormat = null;
//TODO: remove !!!
globalThis.addEventListener('languagechange', function (event) {
  cacheLocale = '';
});
var getNumberFormat = function () {
  var locale = Expression.locale;
  if (cacheLocale !== locale) {
    cacheLocale = locale;
    cachedNumberFormat = new Intl.NumberFormat(locale, {maximumFractionDigits: 1/0, minimumFractionDigits: 0, useGrouping: false});
  }
  return cachedNumberFormat;
};

Expression.Integer.prototype.toMathML = function (options) {
  //var d = Expression.toDecimalString(this, options);
  //if (d != undefined) {
  //  return d;
  //}
  var x = this;
  var sign = x.compareTo(Expression.ZERO) < 0 ? '-' : '';
  var abs = x.compareTo(Expression.ZERO) < 0 ? x.negate() : x;
  var s = abs.value.toString();
  var tmp = getNumberFormat().format(s);
  return sign === "-" ? "<mrow>" + "<mo>&minus;</mo>" + "<mn>" + tmp + "</mn>" + "</mrow>" : "<mn>" + tmp + "</mn>";
};
Expression.BinaryOperation.prototype.toMathML = function (options) {
  options = options == null ? {} : options;
  if (options != null &&
      options.rounding != null &&
      this.unwrap() instanceof Expression.Exponentiation &&
      this.unwrap().a.unwrap() instanceof Expression.Symbol &&
      this.unwrap().a.unwrap() !== Expression.E &&
      this.unwrap().a.unwrap() !== Expression.PI &&
      (this.unwrap().b.unwrap() instanceof Expression.Integer || this.unwrap().b.unwrap() instanceof Expression.Negation && this.unwrap().b.unwrap().b.unwrap() instanceof Expression.Integer)) {
    options = Object.assign({}, options, {rounding: null});
  }
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }

  //!2019-05-16
  if (this instanceof Expression.Addition && options.printId == undefined) {
    var s = [];
    var b = null;
    for (var x of this.summands()) {
      if (b != null) {
        var n = false;
        if (b.isNegative()) {
          n = true;
          b = b.negateCarefully();
        }
        var fence = this.getPrecedence() >= getPrecedence(b, options);
        fence = fence || b.isUnaryPlusMinus();
        s.push((fence ? '<mrow><mo>(</mo>' : '') + b.toMathML(Expression.setTopLevel(fence, options)) + (fence ? '<mo>)</mo></mrow>' : ''));
        s.push(n ? '<mo>&minus;</mo>' : '<mo>+</mo>');
      }
      b = x;
    }
    s = s.reverse().join('');
    var a = b;
    var fence = a.getPrecedence() + (a.isRightToLeftAssociative() ? -1 : 0) < this.getPrecedence();
    if (options != undefined && options.isTopLevel != undefined && options.isTopLevel === false) {
      fence = fence || a.isUnaryPlusMinus();
    }
    s = (fence ? "<mrow><mo>(</mo>" : "") + a.toMathML(Expression.setTopLevel(fence || options == undefined || options.isTopLevel, options)) + (fence ? "<mo>)</mo></mrow>" : "") + s;
    return '<mrow>' + s + '</mrow>';
  }

  //!
  var a = this.a;
  var b = this.b;
  var isSubtraction = false;
  // TODO: check
  if (this instanceof Expression.Addition && b.isNegative()) {
    isSubtraction = true;
    b = b.negateCarefully();//?
  }

  var fa = a.getPrecedence() + (a.isRightToLeftAssociative() ? -1 : 0) < this.getPrecedence();
  var fb = this.getPrecedence() + (this.isRightToLeftAssociative() ? -1 : 0) >= b.getPrecedence();
  if (options != undefined && options.isTopLevel != undefined && options.isTopLevel === false) {
    fa = fa || a.isUnaryPlusMinus();
  }
  fb = fb || b.isUnaryPlusMinus();
  fb = fb || (this.unwrap() instanceof Expression.Exponentiation && b.unwrap() instanceof Expression.Exponentiation);// 2^3^4
  fa = fa || (this.unwrap() instanceof Expression.Exponentiation && a.unwrap() instanceof Expression.Function); // cos(x)^(2+3)
  var s = isSubtraction ? "-" : this.getS();

  if (this instanceof Expression.Exponentiation) {
    if (a.unwrap() === Expression.E && b.unwrap() instanceof Expression.Matrix) {
      return '<mrow><mi>exp</mi><mo>&af;</mo>' + b.toMathML(options) + '</mrow>';
    }
    var boptions = options;
    if (!(a.unwrap() instanceof Expression.Matrix)) {
      boptions = Object.assign({}, options || {}, {nofractions: true});
    }

      return "<msup>" +
             (fa ? "<mrow><mo>(</mo>" : "") + a.toMathML(Expression.setTopLevel(fa || options == undefined || options.isTopLevel, options)) + (fa ? "<mo>)</mo></mrow>" : "") +
             (fb ? "<mrow><mo>(</mo>" : "") + b.toMathML(Expression.setTopLevel(fb, boptions)) + (fb ? "<mo>)</mo></mrow>" : "") +
             "</msup>";
  }
  if (this.isNegation()) {
    // assert(fa === false);
      return "<mrow><mo>&minus;</mo>" + (fb ? "<mrow><mo>(</mo>" : "") + b.toMathML(Expression.setTopLevel(fb, options)) + (fb ? "<mo>)</mo></mrow>" : "") + "</mrow>";
  }
  //TODO: fix spaces (matrix parsing)
  // &times; looks better than &it; when multiplying matrices (?)
  // &sdot; looks better than &times;
  var isScalarOrMatrixSymbol = function (e) {
    return Expression.isScalar(e.unwrap()) || e.unwrap() instanceof Expression.MatrixSymbol;
  };
  var base = function (e) {
    return e instanceof Expression.Exponentiation && e.b.unwrap() instanceof Expression.Integer ? e.a.unwrap() : e;
  };
  var canUseInvisibleTimes = function (e) {
    return !fa && !fb && (e.a.unwrap() instanceof Expression.Integer || base(e.a.unwrap()) instanceof Expression.Symbol || e.a.unwrap() instanceof Expression.Multiplication && base(e.a.unwrap().b.unwrap()) instanceof Expression.Symbol) && (base(e.b.unwrap()) instanceof Expression.Symbol || options.rounding == null && e.b.unwrap() instanceof Expression.SquareRoot);
    //return options.rounding == null && !fa && !fb && isScalarOrMatrixSymbol(e.a) && isScalarOrMatrixSymbol(e.b) && !(e.a instanceof Expression.Integer && (e.b instanceof Expression.Integer || e.b instanceof Expression.Exponentiation && e.b.a instanceof Expression.Integer))
  };
  //!2020-08-02
  if (this instanceof Expression.Multiplication && options.printId == undefined) {
    var f = true;
    for (var x = this; x != null; x = x instanceof Expression.Multiplication ? x.a.unwrap() : null) {
      var factor = x instanceof Expression.Multiplication ? x.b.unwrap() : x;
      if (!(base(factor).unwrap() instanceof Expression.Symbol && Expression.isScalar(factor))) {
        f = false;
      }
    }
    if (f) {
      var s = [];
      for (var x = this; x != null; x = x instanceof Expression.Multiplication ? x.a.unwrap() : null) {
        var factor = x instanceof Expression.Multiplication ? x.b.unwrap() : x;
        s.push(factor.toMathML(Expression.setTopLevel(fence, options)));
      }
      s = s.reverse().join('<mo>&it;</mo>');
      return '<mrow>' + s + '</mrow>';
    }
  }
  return "<mrow>" +
         (fa ? "<mrow><mo>(</mo>" : "") + a.toMathML(Expression.setTopLevel(fa || options == undefined || options.isTopLevel, options)) + (fa ? "<mo>)</mo></mrow>" : "") +
         (s === '*' ? (canUseInvisibleTimes(this) ? '<mo>&it;</mo>' : '<mo>&sdot;</mo>') : (s === '-' ? '<mo>&minus;</mo>' : (s === '/' ? '<mo>&#x2215;</mo>' : (this instanceof Expression.Comma ? '<mo lspace="0em" rspace="0.55em">' + ',' + '</mo>' : '<mo>' + s + '</mo>')))) +
         (fb ? "<mrow><mo>(</mo>" : "") + b.toMathML(Expression.setTopLevel(fb, options)) + (fb ? "<mo>)</mo></mrow>" : "") +
         "</mrow>";
};

Expression.Symbol.prototype.toMathML = function (options) {
  var x = this;
  var s = x.symbol;
  var i = s.indexOf("_");
  if (i !== -1) {
    var indexes = s.slice(i + 1).replace(/^\(|\)$/g, "").split(",");
    var indexesMathML = "";
    for (var j = 0; j < indexes.length; j += 1) {
      indexesMathML += j !== 0 ? "<mo>,</mo>" : "";
      indexesMathML += /^\d+$/.exec(indexes[j]) != undefined ? Expression.Integer.fromString(indexes[j]).toMathML() : "<mi>" + indexes[j] + "</mi>";
    }
    if (indexes.length > 1) {
      indexesMathML = "<mrow>" + indexesMathML + "</mrow>";
    }
    return "<msub>" +
           "<mi>" + s.slice(0, i) + "</mi>" +
           indexesMathML +
           "</msub>";
  }
  return "<mi>" + (this instanceof Expression.IdentityMatrix ? '<span class="dotted-underline" title="' + i18n.identityMatrix + '" aria-label="' + i18n.identityMatrix + '">' + s + '</span>' : s) + "</mi>";
};
Expression.Negation.prototype.toMathML = function (options) {
  var b = this.b;
  var fb = this.getPrecedence() + (this.isRightToLeftAssociative() ? -1 : 0) >= b.getPrecedence();
  fb = fb || b.isUnaryPlusMinus();
  // assert(fa === false);
  return "<mrow><mo>&minus;</mo>" + (fb ? "<mrow><mo>(</mo>" : "") + b.toMathML(Expression.setTopLevel(fb, options)) + (fb ? "<mo>)</mo></mrow>" : "") + "</mrow>";
};

Expression.Factorial.prototype.toMathML = function (options) {
  var n = this.n.unwrap();
  var fn = !(n instanceof Expression.Integer && n.compareTo(Expression.ZERO) > 0);
  return "<mrow>" + (fn ? "<mrow><mo>(</mo>" : "") + n.toMathML(Expression.setTopLevel(fn, options)) + (fn ? "<mo>)</mo></mrow>" : "") + "<mo>!</mo></mrow>";
};

Condition.prototype.toMathML = function (options) {
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
    var c = this.array[i];
    s += (i !== 0 ? '<mo lspace="0em" rspace="0.55em">,</mo>' : '');
    s += '<mrow>';
    var operator = '<mo>' + (c.operator === Condition.NEZ ? '&ne;' : (c.operator === Condition.EQZ ? '=' : (c.operator === Condition.GTZ ? '&gt;' : '???'))) + '</mo>';
    if (c.expression instanceof Expression.Addition && c.expression.a instanceof Expression.Symbol && c.expression.b instanceof Expression.Integer) {//TODO: ?
      var left = c.expression.a;
      var right = c.expression.b.negate();
      s += left.toMathML(options) + operator + right.toMathML(options);
    } else {
      s += c.expression.toMathML(options) + operator + Expression.ZERO.toMathML();
    }
    s += '</mrow>';
  }
  return this.array.length === 1 ? s : '<mrow>' + s + '</mrow>';
};

Expression.Complex.prototype.toMathML = function (options) {
  return this.toStringInternal(options, "<mo>&it;</mo>", "<mi>&ii;</mi>", "<mo>&minus;</mo>", "<mo>+</mo>", "<mrow>", "</mrow>", function (x, options) { return x.toMathML(options); });
};

Expression.GF2.prototype.toMathML = function (options) {
  //TODO: fix
  return this.a.toMathML(options);
};
Expression.GF2Value.prototype.toMathML = function (options) {
  return "<mrow>" + "<mn>" + this.value.toString() + "</mn>" + "</mrow>";
};
Expression.Degrees.prototype.toMathML = function (options) {
  return "<mrow>" + this.value.toMathML(options) + "<mo>&it;</mo><mi>&deg;</mi></mrow>";
};
Expression.Radians.prototype.toMathML = function (options) {
  return "<mrow>" + this.value.toMathML(options) + "<mo>&it;</mo><mi>rad</mi></mrow>";
};

NonSimplifiedExpression.prototype.toMathML = function (options) {
  //?
  //options = options.rounding != null ? Object.assign({}, options, {rounding: null}) : options;
  if (options != null && options.printId != undefined) {
    return "<mrow id=\"" + this.getId() + "\">" + this.e.toMathML(options) + "</mrow>";
  }
  return this.e.toMathML(options);
};

Expression.prototype.toMathML = function (options) {
  throw new TypeError();
};
