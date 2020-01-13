/*jslint plusplus: true, vars: true, indent: 2 */
  import Expression from './Expression.js';
  import Matrix from './Matrix.js';

  //var isAlpha = function (code) {
  //  return (code >= "a".charCodeAt(0) && code <= "z".charCodeAt(0)) ||
  //         (code >= "A".charCodeAt(0) && code <= "Z".charCodeAt(0));
  //};

  // http://en.wikipedia.org/wiki/Operators_in_C_and_C%2B%2B#Operator_precedence

  var LEFT_TO_RIGHT = 0;
  var RIGHT_TO_LEFT = 1;

  var EQUALITY_PRECEDENCE = 1;
  var ADDITIVE_PRECEDENCE = 2;
  var MULTIPLICATIVE_PRECEDENCE = 3;
  var UNARY_PRECEDENCE = 5;

  var UNARY_PRECEDENCE_PLUS_ONE = UNARY_PRECEDENCE + 1; // TODO: remove

  var Operator = function (name, arity, rightToLeftAssociative, precedence, i) {
    this.name = name;
    this.arity = arity;
    this.rightToLeftAssociative = rightToLeftAssociative;
    this.precedence = precedence;
    this.i = i;
    //this.xyz = isAlpha(name.charCodeAt(0)) && isAlpha(name.charCodeAt(name.length - 1));
  };

  var ADDITION = new Operator("+", 2, LEFT_TO_RIGHT, ADDITIVE_PRECEDENCE, function (a, b) {
    return a.add(b);
  });
  var MULTIPLICATION = new Operator("*", 2, LEFT_TO_RIGHT, MULTIPLICATIVE_PRECEDENCE, function (a, b) {
    return a.multiply(b);
  });
  // Exponentiation has precedence as unary operators
  var EXPONENTIATION = new Operator("^", 2, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a, b) {
    return a.pow(b);
  });

  var UNARY_PLUS = new Operator("+", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (e) {
    return e;
  });
  var UNARY_MINUS = new Operator("-", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (e) {
    return e.negate();
  });

  var prepareTrigonometricArgument = function (a) {
    if (a instanceof Expression.Integer) {
      return new Expression.Degrees(a);
    }
    if (a instanceof Expression.NonSimplifiedExpression && a.e instanceof Expression.Integer) {
      return new Expression.NonSimplifiedExpression(new Expression.Degrees(a));
    }
    if (a instanceof Expression.NonSimplifiedExpression && a.e instanceof Expression.Negation && a.e.b instanceof Expression.NonSimplifiedExpression && a.e.b.e instanceof Expression.Integer) {
      return new Expression.NonSimplifiedExpression(new Expression.Degrees(a));
    }
    return a;
  };

  var toDegrees = function (a) {
    return a instanceof Expression.NonSimplifiedExpression ? new Expression.NonSimplifiedExpression(new Expression.Degrees(a)) : new Expression.Degrees(a);
  };

  var toRadians = function (a) {
    return a instanceof Expression.NonSimplifiedExpression ? new Expression.NonSimplifiedExpression(new Expression.Radians(a)) : new Expression.Radians(a);
  };

  var notSupported = function (a) {
    throw new TypeError();
  };

  var operations = [
    new Operator("=", 2, LEFT_TO_RIGHT, EQUALITY_PRECEDENCE, function (a, b) {
      return a.transformEquality(b);
    }),
    new Operator(";", 2, LEFT_TO_RIGHT, EQUALITY_PRECEDENCE, function (a, b) {
      throw new RangeError("NotSupportedError");
      //return a.transformStatement(b);
    }),

    ADDITION,
    new Operator("-", 2, LEFT_TO_RIGHT, ADDITIVE_PRECEDENCE, function (a, b) {
      return a.subtract(b);
    }),
    MULTIPLICATION,
    new Operator("/", 2, LEFT_TO_RIGHT, MULTIPLICATIVE_PRECEDENCE, function (a, b) {
      return a.divide(b);
    }),
    new Operator("\\", 2, LEFT_TO_RIGHT, MULTIPLICATIVE_PRECEDENCE, function (a, b) {
      return a.inverse().multiply(b);
    }),
    //new Operator("%", 2, LEFT_TO_RIGHT, MULTIPLICATIVE_PRECEDENCE, function (a, b) {
    //  return a.remainder(b);
    //}),
    //UNARY_PLUS,
    //UNARY_MINUS,
    EXPONENTIATION,
    new Operator("**", EXPONENTIATION.arity, EXPONENTIATION.rightToLeftAssociative, EXPONENTIATION.precedence, EXPONENTIATION.i),
    new Operator(".^", 2, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a, b) {
      return a.elementWisePower(b);
    }),//?
    new Operator("\u221A", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a) {
      return a.squareRoot();
    }),
    new Operator("sqrt", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      return a.squareRoot();
    }),
    new Operator("\u221B", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a) {
      return a.cubeRoot();
    }),
    new Operator("cbrt", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      return a.cubeRoot();
    }),
    new Operator("\u221C", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a) {
      return a.pow(Expression.ONE.divide(Expression.TWO.add(Expression.TWO)));
    }),
    new Operator("rank", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a) {
      return a.rank();
    }),
    new Operator("adjugate", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a) {
      return a.adjugate();
    }),
    //new Operator("trace", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a) {
    //  return Expression.transformTrace(a);
    //}),
    new Operator("inverse", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a) {
      return a.inverse();
    }),
    new Operator("det", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a) {//?
      return a.determinant();
    }),
    new Operator("determinant", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a) {
      return a.determinant();
    }),
    new Operator("row-reduce", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a) {
      return a.rowReduce();
    }),
    new Operator("transpose", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a) {
      return a.transpose();
    }),
    //new Operator("^T", 1, LEFT_TO_RIGHT, UNARY_PRECEDENCE, function (a) {
    //  return a.transpose();
    //}),
    //new Operator("^t", 1, LEFT_TO_RIGHT, UNARY_PRECEDENCE, function (a) {
    //  return a.transpose();
    //}),
    new Operator("'", 1, LEFT_TO_RIGHT, UNARY_PRECEDENCE, function (a) {
      return a.transpose();
    }),

    //?
    new Operator("solve", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a) {
      return a.transformNoAnswerExpression("solve");//?
    }),

    new Operator("GF2", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a) {
      return a.GF2();
    }),

    new Operator("cos", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      a = prepareTrigonometricArgument(a);
      return a.cos();
    }),
    new Operator("sin", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      a = prepareTrigonometricArgument(a);
      return a.sin();
    }),
    new Operator("tan", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      //a = prepareTrigonometricArgument(a);
      //return a.sin().divide(a.cos());
      var a2 = prepareTrigonometricArgument(a.multiply(Expression.TWO));
      return a2.sin().divide(a2.cos().add(Expression.ONE));
    }),
    new Operator("cot", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      if (a instanceof Expression.Matrix) {
        a = prepareTrigonometricArgument(a);
        return a.cos().divide(a.sin());
      }
      //a = prepareTrigonometricArgument(a);
      //return a.cos().divide(a.sin());
      var a2 = prepareTrigonometricArgument(a.multiply(Expression.TWO));
      return a2.cos().add(Expression.ONE).divide(a2.sin());
    }),
    new Operator("°", 1, LEFT_TO_RIGHT, UNARY_PRECEDENCE, function (a) {
      var x = toDegrees(a);
      if (x == null) {
        throw new RangeError("NotSupportedError");
      }
      return x;
    }),
    new Operator("rad", 1, LEFT_TO_RIGHT, UNARY_PRECEDENCE, function (a) {
      var x = toRadians(a);
      if (x == null) {
        throw new RangeError("NotSupportedError");
      }
      return x;
    }),

    new Operator("exp", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      return a.exp();
    }),
    new Operator("log", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      throw new RangeError("NotSupportedError");
    }),

    new Operator("\\left", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      return a;
    }),
    new Operator("\\right", 1, LEFT_TO_RIGHT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      return a;
    }),

    new Operator("cosh", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      return a.exp().add(a.negate().exp()).divide(Expression.TWO);
    }),
    new Operator("sinh", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      return a.exp().subtract(a.negate().exp()).divide(Expression.TWO);
    }),
    new Operator("tanh", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      return a.exp().subtract(a.negate().exp()).divide(a.exp().add(a.negate().exp()));
    }),

    new Operator("arccos", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, notSupported),
    new Operator("arcsin", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, notSupported),
    new Operator("arctan", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, notSupported),
    new Operator("arcosh", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, notSupported),
    new Operator("arsinh", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, notSupported),
    new Operator("artanh", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, notSupported),

    new Operator("frac", 2, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a, b) {
      return a.divide(b);
    }),

    new Operator("!", 1, LEFT_TO_RIGHT, UNARY_PRECEDENCE, function (a) {
      a = a.unwrap();
      if (!(a instanceof Expression.Integer)) {
        throw new TypeError();
      }
      //return a.factorial();
      if (a.compareTo(Expression.ZERO) < 0) {
        throw new TypeError();
      }
      var f = Expression.ONE;
      for (var i = a; i.compareTo(Expression.ONE) >= 0; i = i.subtract(Expression.ONE)) {
        f = f.multiply(i);
      }
      return f;
    }),
    new Operator("!!", 1, LEFT_TO_RIGHT, UNARY_PRECEDENCE, notSupported), // to not parse 3!! as (3!)!, see https://en.wikipedia.org/wiki/Double_factorial
    new Operator("!!!", 1, LEFT_TO_RIGHT, UNARY_PRECEDENCE, notSupported)
  ];

  function OperationSearchCache() {
    this.map = {};
    this.re = null;
  }

  OperationSearchCache.prototype.append = function (operator) {
    this.map[operator.name.toLowerCase()] = operator;
    this.re = null;
  };
  OperationSearchCache.prototype.getByName = function (name) {
    return this.map[name.toLowerCase()];
  };
  OperationSearchCache.prototype.getRegExp = function () {
    // https://stackoverflow.com/questions/3561493/is-there-a-regexp-escape-function-in-javascript
    var escapeRegExp = function (s) {
      return s.replace(/[-\/\\^$*+?.()|[\]{}]/g, '\\$&');
    };
    // longest: ? "^T" and "^"
    // ignore case
    if (this.re == null) {//TODO: ?
      var names = [];
      for (var name in this.map) {
        if (Object.prototype.hasOwnProperty.call(this.map, name)) {
          names.push(name);
        }
      }
      names.sort(function (a, b) {
        return a + '\uFFFF' < b + '\uFFFF' ? -1 : +1;
      });
      var source = new Array(names.length);
      for (var i = 0; i < names.length; i += 1) {
        source[i] = escapeRegExp(names[i]);
      }
      this.re = new RegExp('^(?:' + source.join('|') + ')', 'i');
    }
    return this.re;
  };

  var operationSearchCache = new OperationSearchCache();

  var i = -1;
  while (++i < operations.length) {
    operationSearchCache.append(operations[i]);
  }

  var nextToken = function (tokenizer) {
    var token = null;
    do {
      token = tokenizer.next();
    } while (token.type === 'whitespace');
    return token;
  };

  var parsePunctuator = function (tokenizer, token, punctuator) {
    if (token.type !== 'punctuator' || token.value !== punctuator) {
      ExpressionParser.startPosition = tokenizer.position - token.value.length;
      ExpressionParser.endPosition = tokenizer.position;
      ExpressionParser.input = tokenizer.input;
      if (token.type === 'EOF') {
        throw new RangeError("UserError: unexpected end of input, '" + punctuator + "' expected");
      }
      throw new RangeError("UserError: unexpected '" + token.value + "', '" + punctuator + "' expected");
    }
    token = nextToken(tokenizer);
    return token;
  };

  function ParseResult(result, token) {
    this.result = result;
    this.token = token;
  }

  var parseMatrix = function (tokenizer, token, context) {
    var openingBracket = "{";
    var closingBracket = "}";

    var rows = [];
    var hasNextRow = true;
    while (hasNextRow) {
      token = parsePunctuator(tokenizer, token, openingBracket);
      var row = [];
      var hasNextCell = true;
      while (hasNextCell) {
        var tmp = parseExpression(tokenizer, token, context, 0, undefined);
        token = tmp.token;
        row.push(tmp.result);
        if (token.type === 'punctuator' && token.value === ",") {
          hasNextCell = true;
          token = nextToken(tokenizer);
        } else {
          hasNextCell = false;
        }
      }
      token = parsePunctuator(tokenizer, token, closingBracket);
      rows.push(row);
      if (token.type === 'punctuator' && token.value === ",") {
        hasNextRow = true;
        token = nextToken(tokenizer);
      } else {
        hasNextRow = false;
      }
    }
    token = parsePunctuator(tokenizer, token, "}");
    return new ParseResult(context.wrap(Expression.Matrix.fromArray(rows)), token);
  };

  var parseLaTeXMatrix = function (tokenizer, token, context, rowSeparator) {
    var rows = [];
    var firstRow = true;
    while (firstRow || (token.type === 'punctuator' && token.value === rowSeparator)) {
      if (firstRow) {
        firstRow = false;
      } else {
        token = nextToken(tokenizer);
      }
      var row = [];
      var firstCell = true;
      while (firstCell || token.type === 'punctuator' && token.value === "&") {
        if (firstCell) {
          firstCell = false;
        } else {
          token = nextToken(tokenizer);
        }
        var tmp = parseExpression(tokenizer, token, context, ADDITIVE_PRECEDENCE - 1, undefined);
        token = tmp.token;
        row.push(tmp.result);
      }
      rows.push(row);
    }
    return new ParseResult(context.wrap(Expression.Matrix.fromArray(rows)), token);
  };

  var parseLaTeXArgument = function (tokenizer, token, context) {
    return parseExpression(tokenizer, token, context, 0, undefined);
  };

  var getVulgarFraction = function (vulgarFraction) {
    var input = normalizeVulgarFractions(vulgarFraction);
    var e = Expression.Integer.fromString(input.slice(0, input.indexOf('/'))).divide(Expression.Integer.fromString(input.slice(input.indexOf('/') + '/'.length)));
    return e;
  };

  var getDecimalFraction = function (integerPart, nonRepeatingFractionalPart, repeatingFractionalPart, exponentPart) {
    var numerator = Expression.ZERO;
    var denominator = Expression.ONE;

    if (integerPart != undefined) {
      numerator = Expression.Integer.fromString(integerPart);
    }
    if (nonRepeatingFractionalPart != undefined) {
      var factor = Expression.pow(Expression.TEN, nonRepeatingFractionalPart.length);
      numerator = numerator.multiply(factor).add(Expression.Integer.fromString(nonRepeatingFractionalPart));
      denominator = denominator.multiply(factor);
    }
    if (repeatingFractionalPart != undefined) {
      var factor = Expression.pow(Expression.TEN, repeatingFractionalPart.length).subtract(Expression.ONE);
      numerator = numerator.multiply(factor).add(Expression.Integer.fromString(repeatingFractionalPart));
      denominator = denominator.multiply(factor);
    }
    if (exponentPart != undefined) {
      var exponent = 0 + Number.parseInt(exponentPart, 10);
      var factor = Expression.pow(Expression.TEN, exponent < 0 ? -exponent : exponent);
      if (exponent < 0) {
        denominator = denominator.multiply(factor);
      } else {
        numerator = numerator.multiply(factor);
      }
    }

    var value = numerator.divide(denominator);
    return value;
  };

  var parseDecimalFraction = function (tokenizer, token, context) {
    var isOnlyInteger = true;
    var result = undefined;
    if (token.type === 'integerLiteral') {
      result = Expression.Integer.fromString(token.value);
      result = context.wrap(result);
      token = nextToken(tokenizer);
    } else if (token.type === 'numericLiteral') {
      var value = token.value;
      //var match = token.match;
      var match = decimalFractionWithGroups.exec(value);
      isOnlyInteger = false;
      result = getDecimalFraction(match[1], match[2], match[3], match[4]);
      result = context.wrap(result);
      token = nextToken(tokenizer);
    }
    //!
    if (isOnlyInteger || result == undefined) {
      if (token.type === 'vulgarFraction') {
        var fraction = context.wrap(getVulgarFraction(token.value, context));
        if (result != undefined) {
          result = ADDITION.i(result, fraction).addPosition(tokenizer.position - token.value.length, ADDITION.name.length, tokenizer.input);
        } else {
          result = fraction;
        }
        token = nextToken(tokenizer);
      }
    }
    return result != undefined ? new ParseResult(result, token) : undefined;
  };

  // TODO: sticky flags - /\s+/y
  var whiteSpaces = /^\s+/;
  var punctuators = /^(?:[,&(){}|■@]|\\\\|(?:\\begin|\\end)(?:\{[bvp]?matrix\})?)/;
  var integerLiteral = /^\d+(?![\d.,eE])/; // for performance
  var integerLiteralWithoutComma = /^\d+(?![\d.eE])/; // for performance
  var decimalFraction = /^(?=[.,]?\d)\d*(?:[.,]\d*(?:\(\d+\))?)?(?:[eE][\+\-]?\d+)?/;
  var decimalFractionWithoutComma = /^(?=[.]?\d)\d*(?:[.]\d*(?:\(\d+\))?)?(?:[eE][\+\-]?\d+)?/;
  // Base Latin, Base Latin upper case, Base Cyrillic, Base Cyrillic upper case, Greek alphabet
  var symbols = /^(?:alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|mu|nu|xi|omicron|pi|rho|varsigma|sigma|tau|upsilon|phi|chi|psi|omega|circ|[a-zA-Z\u0430-\u044F\u0410-\u042F\u03B1-\u03C9])(?:\_\d+|\_\([a-z\d]+,[a-z\d]+\)|[\u2080-\u2089]+)?/;
  var superscripts = /^[\u00B2\u00B3\u00B9\u2070\u2074-\u2079]+/; // superscript digits 2310456789
  var vulgarFractions = /^[\u00BC-\u00BE\u2150-\u215E]/;
  //var other = /^\S/u;
  var other = /^(?:[\uD800-\uDBFF][\uDC00-\uDFFF]|\S)/; // should not split surrogate pairts (for Tokenizer and other things)

  var decimalFractionWithGroups = /^(\d+)?(?:[.,](\d+)?(?:\((\d+)\))?)?(?:[eE]([\+\-]?\d+))?$/;

  // s.normalize("NFKD").replace(/[\u2044]/g, "/")
  var normalizeSuperscripts = function (s) {
    return s.replace(/[\u00B2\u00B3\u00B9\u2070\u2074-\u2079]/g, function (c) {
      var charCode = c.charCodeAt(0);
      if (charCode === 0x00B2) {
        return "2";
      }
      if (charCode === 0x00B3) {
        return "3";
      }
      if (charCode === 0x00B9) {
        return "1";
      }
      return (charCode - 0x2074 + 4).toString();
    });
  };

  var normalizeVulgarFractions = function (s) {
    return s.replace(/[\u00BC-\u00BE\u2150-\u215E]/g, function (c) {
      var charCode = c.charCodeAt(0);
      var i = charCode - 0x2150 < 0 ? (charCode - 0x00BC) * 2 : (3 + charCode - 0x2150) * 2;
      return "141234171911132315253545165618385878".slice(i, i + 2).replace(/^\S/g, "$&/").replace(/1\/1/g, "1/10");
    });
  };

  var normalizeSubscripts = function (s) {
    var i = s.length - 1;
    while (i >= 0 && s.charCodeAt(i) >= 0x2080) {
      i -= 1;
    }
    return i === s.length - 1 ? s : s.slice(0, i + 1) + "_" + s.slice(i + 1).replace(/[\u2080-\u2089]/g, function (c) {
      return String.fromCharCode(c.charCodeAt(0) - 0x2080 + "0".charCodeAt(0));
    });
  };

  var normalizeGreek = function (s) {
    var i = s.indexOf("_");
    var k = i === -1 ? s.length : i;
    if (k > 1) {
      var name = s.slice(0, k);
      var greek = " alpha beta gamma delta epsilon zeta eta theta iota kappa lambda mu nu xi omicron pi rho varsigma sigma tau upsilon phi chi psi omega ";
      var j = greek.indexOf(" " + name + " ");
      if (j !== -1) {
        return String.fromCharCode(0x03B1 + greek.slice(0, j).split(" ").length - 1) + s.slice(k);
      }
    }
    return s;
  };

  var parseExpression = function (tokenizer, token, context, precedence, left) {
    var ok = true;
    var isDecimalFraction = false;
    var tmp = undefined;
    var right = undefined;
    //!

    while (token.type !== 'EOF' && ok) {
      var op = undefined;
      var operand = undefined;

        var bestMatch = token.type === 'operator' ? operationSearchCache.getByName(token.value) : null;
        if (bestMatch != null) {
          op = left == null && bestMatch.name === '+' ? UNARY_PLUS : (left == null && bestMatch.name === '-' ? UNARY_MINUS : bestMatch);
        }
        //  if (Input.startsWith(input, position, '\\begin') || Input.startsWith(input, position, '\\end')) {
        //    op = null;
        //  }

      //if (op != null && op.name === "\\" && Input.startsWith(input, position, "\\\\")) {
      //  if (isMatrixElement) {//TODO: optimize
      //    op = null;
        //} else if (Input.startsWith(input, position + 1, "begin") || Input.startsWith(input, position + 1, "left")) {
        //  op = null;
      //  }
      //}

      if (op != null && op.name === "frac") { // !isAlpha(Input.getFirst(input, position + "frac".length))
        if (!(left == null && precedence <= UNARY_PRECEDENCE_PLUS_ONE || precedence < MULTIPLICATION.precedence)) {
          ok = false;
        } else {
        // https://en.wikipedia.org/wiki/Operand#Positioning_of_operands - prefix notation

        token = nextToken(tokenizer);
        tmp = parseExpression(tokenizer, token, context, MULTIPLICATION.precedence, undefined);
        var a = tmp.result;
        token = tmp.token;
        tmp = parseExpression(tokenizer, token, context, MULTIPLICATION.precedence, undefined);
        var b = tmp.result;
        token = tmp.token;
        // addPosition - ?
        operand = op.i(a, b);
        ok = true;
        }
      } else if (op != undefined) {
        // TODO: check if the checks are needed (tests - ?)
        if (!(left != undefined && (op.arity !== 1 || op.rightToLeftAssociative !== RIGHT_TO_LEFT || precedence < MULTIPLICATION.precedence) ||
              left == undefined && op.arity === 1 && op.rightToLeftAssociative === RIGHT_TO_LEFT) ||
            //!(!candidate.xyz || !isAlpha(Input.getFirst(input, position + candidate.name.length))) ||//TODO: fix - ExpressionParser.parse("George")
            precedence > op.precedence + (op.rightToLeftAssociative === RIGHT_TO_LEFT ? 0 : -1)) {
          ok = false;
        } else {
          var operatorPosition = tokenizer.position - token.value.length;
          token = nextToken(tokenizer);
          if (op.arity === 1 && op.rightToLeftAssociative !== RIGHT_TO_LEFT) {
            //TODO: fix
            ExpressionParser.startPosition = operatorPosition;
            ExpressionParser.endPosition = operatorPosition + op.name.length;
            ExpressionParser.input = tokenizer.input;
            left = op.i(left).addPosition(operatorPosition, op.name.length, tokenizer.input);
          } else {
            if (op.arity === 1 && op.rightToLeftAssociative === RIGHT_TO_LEFT && op.precedence === UNARY_PRECEDENCE_PLUS_ONE && op.name.length > 1 &&
                (op.name === "sin" ||
                 op.name === "cos" ||
                 op.name === "sen" ||
                 op.name === "tan" ||
                 op.name === "tg" ||
                 op.name === "cot" ||
                 op.name === "ctg") &&
                (token.type === 'operator' && token.value === EXPONENTIATION.name || token.type === 'superscript')) {
              // https://en.wikipedia.org/wiki/Exponentiation#Exponential_notation_for_function_names

              // cos^2(x)
              //!new 2017-11-04
              // parse an operator for the exponentiation
              var exponentiationPosition = tokenizer.position;

              var exponentiationLength = 0;
              var middle = null;
              if (token.type === 'superscript') {
                var superscript = token.value;
                exponentiationLength = token.value.length;
                token = nextToken(tokenizer);
                middle = Expression.Integer.fromString(normalizeSuperscripts(superscript));
              } else {
                exponentiationLength = EXPONENTIATION.name.length;
                token = nextToken(tokenizer);
                if (token.type !== 'integerLiteral') {
                  ok = false;
                } else {
                  tmp = parseExpression(tokenizer, token, context, EXPONENTIATION.precedence, undefined);
                  middle = tmp.result;
                  token = tmp.token;
                }
              }
              if (ok) {
                // parse an operator for the current operator
                tmp = parseExpression(tokenizer, token, context, op.precedence, undefined);
                right = tmp.result;
                token = tmp.token;
                operand = EXPONENTIATION.i(op.i(right).addPosition(operatorPosition, op.name.length, tokenizer.input), middle).addPosition(exponentiationPosition, exponentiationLength, tokenizer.input);
              }
            } else {
              tmp = parseExpression(tokenizer, token, context, op.precedence, undefined);
              right = tmp.result;
              token = tmp.token;
              //TODO: fix `1/(2-2)`
              ExpressionParser.startPosition = operatorPosition;
              ExpressionParser.endPosition = operatorPosition + op.name.length;
              ExpressionParser.input = tokenizer.input;
              if (op.arity === 1) {
                // left <implicit multiplication> operand
                operand = op.i(right).addPosition(operatorPosition, op.name.length, tokenizer.input);
              } else if (op.arity === 2) {
                left = op.i(left, right).addPosition(operatorPosition, op.name.length, tokenizer.input);
              } else {
                throw new RangeError();
              }
            }
          }
        }
      } else if (left == undefined || precedence < MULTIPLICATION.precedence || (precedence === UNARY_PRECEDENCE_PLUS_ONE && isDecimalFraction && token.type === 'symbol')) {
        if ((tmp = parseDecimalFraction(tokenizer, token, context)) != undefined) {
          operand = tmp.result;
          token = tmp.token;
          isDecimalFraction = true;
        } else if (token.type === 'punctuator' && token.value === "(") {
          token = parsePunctuator(tokenizer, token, "(");
          tmp = parseExpression(tokenizer, token, context, 0, undefined);
          operand = tmp.result;
          token = tmp.token;
          token = parsePunctuator(tokenizer, token, ")");
        } else if (token.type === 'punctuator' && token.value === "{") {
          token = parsePunctuator(tokenizer, token, "{");
          if (token.type === 'punctuator' && token.value === "{") {
            tmp = parseMatrix(tokenizer, token, context);
            operand = tmp.result;
            token = tmp.token;
          } else {
            tmp = parseLaTeXArgument(tokenizer, token, context);
            operand = tmp.result;
            token = tmp.token;
            token = parsePunctuator(tokenizer, token, "}");
          }
        } else if (token.type === 'punctuator' && (token.value === "\\begin{bmatrix}" ||
                                                   token.value === "\\begin{vmatrix}" ||
                                                   token.value === "\\begin{pmatrix}" ||
                                                   token.value === "\\begin{matrix}")) {
          var kind = token.value.slice('\\begin{'.length, -1);
          token = nextToken(tokenizer);
          tmp = parseLaTeXMatrix(tokenizer, token, context, '\\\\');
          operand = tmp.result;
          token = tmp.token;
          if (token.type === 'punctuator' && token.value === "\\end{" + kind + "}") {
            token = nextToken(tokenizer);
          }
          if (kind === 'vmatrix') {
            operand = operand.determinant();//!
          }
        } else if (token.type === 'symbol') {
          var symbolName = token.value;
          symbolName = normalizeSubscripts(symbolName);
          symbolName = normalizeGreek(symbolName);
          operand = context.get(symbolName);
          operand = context.wrap(operand);
          token = nextToken(tokenizer);
        } else if (token.type === 'punctuator' && token.value === "|" && left == undefined) {
          token = parsePunctuator(tokenizer, token, "|");
          tmp = parseExpression(tokenizer, token, context, 0, undefined);
          operand = tmp.result;
          token = tmp.token;
          token = parsePunctuator(tokenizer, token, "|");
          
          operand = operand.determinant();//!
        } else if (token.type === 'punctuator' && token.value === '■') {
          token = parsePunctuator(tokenizer, token, '■');
          token = parsePunctuator(tokenizer, token, '(');
          tmp = parseLaTeXMatrix(tokenizer, token, context, '@');
          operand = tmp.result;
          token = tmp.token;
          token = parsePunctuator(tokenizer, token, ')');
        } else {
          ok = false;
        }
      } else {
        ok = false;
      }

      //!TODO: fix
      if (!ok && left != undefined && precedence <= EXPONENTIATION.precedence + (EXPONENTIATION.rightToLeftAssociative === RIGHT_TO_LEFT ? 0 : -1)) {
        if (token.type === 'superscript') {
          // implicit exponentiation
          //TODO: check position
          var superscript = token.value;
          left = EXPONENTIATION.i(left, Expression.Integer.fromString(normalizeSuperscripts(superscript))).addPosition(tokenizer.position - token.value.length, EXPONENTIATION.name.length, tokenizer.input);
          token = nextToken(tokenizer);
          ok = true;//!
        }
      }

      if (!ok && token.type === 'operator' && token.value === "\\") { // isAlpha(Input.getFirst(input, position + 1))
        // TODO: LaTeX - ?
        ok = true;
        token = nextToken(tokenizer);
      }

      if (operand != undefined) {
        if (left != undefined) {
          // implied multiplication
          var oldPosition = tokenizer.position;
          tmp = parseExpression(tokenizer, token, context, MULTIPLICATION.precedence, operand);
          var right1 = tmp.result;
          token = tmp.token;
          left = MULTIPLICATION.i(left, right1).addPosition(oldPosition, MULTIPLICATION.name.length, tokenizer.input);
        } else {
          left = operand;
        }
      }
    }

    if (left == undefined) {
      ExpressionParser.startPosition = tokenizer.position - token.value.length;
      ExpressionParser.endPosition = tokenizer.position;
      ExpressionParser.input = tokenizer.input;
      if (token.type === 'EOF') {
        throw new RangeError("UserError: unexpected end of input");//TODO: fix
      }
      //TODO: ?
      throw new RangeError("UserError: unexpected '" + token.value + "'");//TODO: fix
    }
    return new ParseResult(left, token);
  };

  var replaceHanidec = function (code) {
    switch (code) {
      case 0x3007: return 0;
      case 0x4E00: return 1;
      case 0x4E8C: return 2;
      case 0x4E09: return 3;
      case 0x56DB: return 4;
      case 0x4E94: return 5;
      case 0x516D: return 6;
      case 0x4E03: return 7;
      case 0x516B: return 8;
      case 0x4E5D: return 9;
    }
    return -1;
  };

  // https://www.ecma-international.org/ecma-402/5.0/index.html#table-numbering-system-digits
  // TODO: remove or add tests
  var replaceSimpleDigit = function (code) {
    if (code >= 0x0660 && code <= 0x0669) {
      return {offset: 0x0660, name: "arab"};
    }
    if (code >= 0x06F0 && code <= 0x06F9) {
      return {offset: 0x06F0, name: "arabext"};
    }
    if (code >= 0x0966 && code <= 0x096F) {
      return {offset: 0x0966, name: "deva"};
    }
    if (code >= 0x09E6 && code <= 0x09EF) {
      return {offset: 0x09E6, name: "beng"};
    }
    if (code >= 0x0A66 && code <= 0x0A6F) {
      return {offset: 0x0A66, name: "guru"};
    }
    if (code >= 0x0AE6 && code <= 0x0AEF) {
      return {offset: 0x0AE6, name: "gujr"};
    }
    if (code >= 0x0B66 && code <= 0x0B6F) {
      return {offset: 0x0B66, name: "orya"};
    }
    if (code >= 0x0BE6 && code <= 0x0BEF) {
      return {offset: 0x0BE6, name: "tamldec"};
    }
    if (code >= 0x0C66 && code <= 0x0C6F) {
      return {offset: 0x0C66, name: "telu"};
    }
    if (code >= 0x0CE6 && code <= 0x0CEF) {
      return {offset: 0x0CE6, name: "knda"};
    }
    if (code >= 0x0D66 && code <= 0x0D6F) {
      return {offset: 0x0D66, name: "mlym"};
    }
    if (code >= 0x0E50 && code <= 0x0E59) {
      return {offset: 0x0E50, name: "thai"};
    }
    if (code >= 0x0ED0 && code <= 0x0ED9) {
      return {offset: 0x0ED0, name: "laoo"};
    }
    if (code >= 0x0F20 && code <= 0x0F29) {
      return {offset: 0x0F20, name: "tibt"};
    }
    if (code >= 0x1040 && code <= 0x1049) {
      return {offset: 0x1040, name: "mymr"};
    }
    if (code >= 0x17E0 && code <= 0x17E9) {
      return {offset: 0x17E0, name: "khmr"};
    }
    if (code >= 0x1810 && code <= 0x1819) {
      return {offset: 0x1810, name: "mong"};
    }
    if (code >= 0x1946 && code <= 0x194F) {
      return {offset: 0x1946, name: "limb"};
    }
    if (code >= 0x1B50 && code <= 0x1B59) {
      return {offset: 0x1B50, name: "bali"};
    }
    if (code >= 0xFF10 && code <= 0xFF19) {
      return {offset: 0xFF10, name: "fullwide"};
    }
    var digit = replaceHanidec(code);
    if (digit !== -1) {
      return {offset: code - digit, name: "hanidec"};
    }
    return undefined;
  };

  var getCharCodeReplacement = function (charCode) {
    if (charCode === 0x0410 || charCode === 0x0430) {
      return "A";
    }
    if (charCode === 0x0412 || charCode === 0x0432) {
      return "B";
    }
    if (charCode === 0x0421 || charCode === 0x0441) {
      return "C";
    }
    //if (charCode === 0x0425 || charCode === 0x0445) {
    //  return "X";
    //}
    if (charCode === 0x0422 || charCode === 0x0442) {
      return "T";
    }
    if (charCode === 0x2212) {
      return "-";
    }
    if (charCode >= 0x2010 && charCode <= 0x2015) {
      return "-";
    }
    if (charCode === 0x00B7 || charCode === 0x00D7 || charCode === 0x2022 || charCode === 0x22C5) {
      return "*";
    }
    // 0x003A - Deutsch
    if (charCode === 0x003A || charCode === 0x00F7) {
      return "/";
    }
    if (charCode >= 0xFF01 && charCode <= 0xFF5E) {
      // normalize full-widht forms:
      return String.fromCharCode(charCode - 0xFF01 + 0x0021);
    }
    if (charCode === 0x060C || charCode === 0x066B) {
      return ",";
    }
    var y = replaceSimpleDigit(charCode);
    if (y != undefined) {
      // TODO: remove
      if (typeof hit === "function") {
        hit({offset: y.name});
      }
      return String.fromCharCode(charCode - y.offset + '0'.charCodeAt(0));
    }
    if (charCode === 0x2147) {
      return "e";
    }
    if (charCode === 0x2148) {
      return "i";
    }
    if (charCode === 0x2215) {
      return "/";
    }
    if (charCode === 0x2713) {
      return "\u221A";
    }
    if (charCode === 0x2061) {
      return " ";
    }
    if (charCode === 0x2062) {
      return "*";
    }
    if (charCode === 0x2063) {
      return ",";
    }
    if (charCode === "〖".charCodeAt(0)) {
      return "(";
    }
    if (charCode === "〗".charCodeAt(0)) {
      return ")";
    }
    if (charCode === "ˆ".charCodeAt(0)) {
      return "^";
    }
    if (charCode === "[".charCodeAt(0)) {//TODO: ?
      return "(";
    }
    if (charCode === "]".charCodeAt(0)) {//TODO: ?
      return ")";
    }
    if (charCode === "∙".charCodeAt(0)) {
      return "*";
    }
    if (charCode === "ー".charCodeAt(0)) {
      return "-";
    }
    return undefined;
  };
  //var replaceRegExp = /[...]/g;
  //var replaceFunction = function (c) {
  //  return getCharCodeReplacement(c.charCodeAt(0));
  //};
  //input = input.replace(replaceRegExp, replaceFunction); - slow in Chrome
  var replaceSomeChars = function (input) {
    var lastIndex = 0;
    var result = '';
    for (var i = 0; i < input.length; i += 1) {
      var charCode = input.charCodeAt(i);
      if (charCode > 0x007F || charCode === 0x003A || charCode === 0x005B || charCode === 0x005D) {
        var x = getCharCodeReplacement(charCode);
        if (x != undefined) {
          result += input.slice(lastIndex, i);
          result += x;
          lastIndex = i + 1;
        }
      }
    }
    result += input.slice(lastIndex);
    return result;
  };

  var config = [
    {type: 'integerLiteral', re: null},
    {type: 'numericLiteral', re: null},
    {type: 'whitespace', re: whiteSpaces},
    {type: 'punctuator', re: punctuators},
    {type: 'operator', re: null},
    {type: 'symbol', re: symbols},
    {type: 'vulgarFraction', re: vulgarFractions},
    {type: 'superscript', re: superscripts},
    {type: 'OTHER', re: other}
  ];

  function Token(type, value) {
    this.type = type;
    this.value = value;
  }

  Token.EOF = new Token('EOF', '', null);

  function Tokenizer(input, position, states) {
    this.input = input;
    this.position = position;
    this.states = states;
  }

  Tokenizer.prototype.next = function () {
    if (this.position >= this.input.length) {
      return Token.EOF;
    }
    // iteration by object keys is slower (?)
    for (var i = 0; i < config.length; i += 1) {
      var c = config[i];
      var type = c.type;
      var re = c.re;
      if (re == null) {
        if (type === 'integerLiteral') {
          if (this.states != null && this.states.value === '{}') {
            re = integerLiteralWithoutComma;
          } else {
            re = integerLiteral;
          }
        } else if (type === 'numericLiteral') {
          if (this.states != null && this.states.value === '{}') {
            re = decimalFractionWithoutComma;
          } else {
            re = decimalFraction;
          }
        } else if (type === 'operator') {
          re = operationSearchCache.getRegExp();//?TODO: 
        }
      }
      var tmp = re.exec(this.input.slice(this.position));
      if (tmp != null) {
        var value = tmp[0];
        if (type === 'punctuator') {
          if (value === '(') {
            this.states = {previous: this.states, value: '()'};
          } else if (value === ')') {
            if (this.states != null && this.states.value === '()') {
              this.states = this.states.previous;
            }
          } else if (value === '{') {
            this.states = {previous: this.states, value: '{}'};
          } else if (value === '}') {
            if (this.states != null && this.states.value === '{}') {
              this.states = this.states.previous;
            }
          }
        }
        this.position += value.length;
        return new Token(type, value);
      }
    }
    throw new TypeError();
  };

  var fs = {};//!TODO: remove!!!

  function ExpressionParser() {
  }

  ExpressionParser.parse = function (input, context) {
    context = context == undefined ? new ExpressionParser.Context(undefined, false) : context;

    ExpressionParser.startPosition = -1;
    ExpressionParser.endPosition = -1;
    ExpressionParser.input = input; //?

    // TODO: remove
    if (typeof input !== "string") {
      throw new RangeError();
    }

    //TODO: fix ???
    input = replaceSomeChars(input);

    if (typeof hit === "function" && context.getter != undefined) {
      var re = /[a-z][a-z][a-z\-]+/gi;
      var m = null;
      while ((m = re.exec(input)) != null) {
        var t = m[0];
        if (!(t in fs) && t.indexOf("-") === -1) {
          fs[t] = true;
          hit({fs: t});
        }
      }
    }

    var tokenizer = new Tokenizer(input, 0, null);
    var token = nextToken(tokenizer);
    var tmp = parseExpression(tokenizer, token, context, 0, undefined);
    token = tmp.token;
    if (token.type !== 'EOF') {
      ExpressionParser.startPosition = tokenizer.position - token.value.length;
      ExpressionParser.endPosition = tokenizer.position;
      ExpressionParser.input = input;
      throw new RangeError("UserError: unexpected '" + token.value + "'");
    }

    return tmp.result;
  };
  
  globalThis.Tokenizer = Tokenizer;

  ExpressionParser.startPosition = -1;
  ExpressionParser.endPosition = -1;
  ExpressionParser.input = "";

  var getConstant = function (symbolName) {
    if (symbolName === "pi" || symbolName === "\u03C0") {
      return Expression.PI;
    }
    if (symbolName === "e") {
      return Expression.E;
    }
    if (symbolName === "i") {
      return Expression.I;
    }
    if (symbolName === "I" || symbolName === "U" || symbolName === "E") {
      return new Expression.IdentityMatrix(symbolName);
    }
    if (symbolName === "circ") { //TODO: ○ - ?
      return Expression.CIRCLE;
    }
    return new Expression.Symbol(symbolName);
  };

  ExpressionParser.Context = function (getter, needsWrap) {
    this.getter = getter;
    this.needsWrap = needsWrap == undefined ? true : needsWrap;
  };
  ExpressionParser.Context.prototype.get = function (symbolName) {
    if (this.getter != undefined) {
      var x = this.getter(symbolName);
      if (x != undefined) {
        return x;
      }
    }
    return getConstant(symbolName);
  };
  ExpressionParser.Context.prototype.wrap = function (e) {
    if (!this.needsWrap) {
      return e;
    }
    return new Expression.NonSimplifiedExpression(e);
  };

  ExpressionParser.addOperation = function (denotation, arity) {
    var newOperation = arity === 1 ? new Operator(denotation, arity, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (a) {
      return a.transformNoAnswerExpression(denotation);
    }) : new Operator(denotation, arity, RIGHT_TO_LEFT, MULTIPLICATIVE_PRECEDENCE, function (a, b) {
      return a.transformNoAnswerExpression(denotation, b);
    });
    //operations.push(newOperation);
    operationSearchCache.append(newOperation);
  };

  ExpressionParser.addDenotations = function (denotationsByOperation) {
    for (var operationName in denotationsByOperation) {
      if (Object.prototype.hasOwnProperty.call(denotationsByOperation, operationName)) {
        var denotations = denotationsByOperation[operationName];
        var operation = operationSearchCache.getByName(operationName);
        var added = {};
        added[operationName] = true;
        for (var key in denotations) {
          if (Object.prototype.hasOwnProperty.call(denotations, key)) {
            var denotation = denotations[key];
            if (added[denotation] == undefined) {
              added[denotation] = true;
              var newOperation = new Operator(denotation, operation.arity, operation.rightToLeftAssociative, operation.precedence, operation.i);
              //operations.push(newOperation);
              operationSearchCache.append(newOperation);
            }
          }
        }
      }
    }
  };

  export default ExpressionParser;
