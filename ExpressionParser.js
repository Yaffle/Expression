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

  var toDegrees = function (a) {
    if (a instanceof Expression.Integer) {
      return new Expression.Degrees(a);
    }
    if (a instanceof Expression.NonSimplifiedExpression && a.e instanceof Expression.Integer) {
      return new Expression.NonSimplifiedExpression(new Expression.Degrees(a));
    }
    if (a instanceof Expression.NonSimplifiedExpression && a.e instanceof Expression.Negation && a.e.b instanceof Expression.NonSimplifiedExpression && a.e.b.e instanceof Expression.Integer) {
      return new Expression.NonSimplifiedExpression(new Expression.Degrees(a));
    }
    return null;
  };

  var prepareTrigonometricArgument = function (a) {
    var x = toDegrees(a);
    return x == null ? a : x;
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
    new Operator("+", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (e) {
      return e;
    }),
    new Operator("-", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE, function (e) {
      return e.negate();
    }),
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

    new Operator("sin", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      if (a.sin == undefined || a.cos == undefined) {
        throw new RangeError("NotSupportedError");
      }
      a = prepareTrigonometricArgument(a);
      return a.sin();
    }),
    new Operator("cos", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      if (a.sin == undefined || a.cos == undefined) {
        throw new RangeError("NotSupportedError");
      }
      a = prepareTrigonometricArgument(a);
      return a.cos();
    }),
    new Operator("tan", 1, RIGHT_TO_LEFT, UNARY_PRECEDENCE_PLUS_ONE, function (a) {
      if (a.sin == undefined || a.cos == undefined) {
        throw new RangeError("NotSupportedError");
      }
      a = prepareTrigonometricArgument(a);
      return a.sin().divide(a.cos());
    }),
    new Operator("°", 1, LEFT_TO_RIGHT, UNARY_PRECEDENCE, function (a) {
      var x = toDegrees(a);
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
    
    new Operator("\\end", 1, RIGHT_TO_LEFT, 0, function (a) {
      throw new RangeError("NotSupportedError");
    })
  ];

  function OperationSearchCache() {
    this.array = new Array(32);
    for (var i = 0; i < this.array.length; i += 1) {
      this.array[i] = undefined;
    }
  }

  var asciiToLowerCase = function (charCode) {
    return (charCode < 65 ? charCode : (charCode < 91 ? charCode + 32 : (charCode < 128 ? charCode : charCode)));
  };

  OperationSearchCache.prototype.append = function (firstCharCode, operator) {
    var index = asciiToLowerCase(firstCharCode) % this.array.length;
    if (this.array[index] == undefined) {
      this.array[index] = [];
    }
    this.array[index].push(operator);
  };
  OperationSearchCache.prototype.getAll = function (firstCharCode) {
    var index = asciiToLowerCase(firstCharCode) % this.array.length;
    return this.array[index];
  };

  var operationSearchCache = new OperationSearchCache();

  var i = -1;
  while (++i < operations.length) {
    operationSearchCache.append(operations[i].name.charCodeAt(0), operations[i]);
  }

  function Input() {
  }

  Input.EOF = -1;
  Input.trimLeft = function (input, position, skip) {
    var match = undefined;
    if ((match = Input.exec(input, position + skip, whiteSpaces)) != undefined) {
      return position + skip + match.length;
    }
    return position + skip;
  };
  Input.parseCharacter = function (input, position, characterCode) {
    var c = Input.getFirst(input, position);
    if (c !== characterCode) {
      ExpressionParser.startPosition = position;
      ExpressionParser.endPosition = position + 1;
      ExpressionParser.input = input;
      if (c === -1) {
        throw new RangeError("UserError: unexpected end of input, '" + String.fromCharCode(characterCode) + "' expected");
      }
      throw new RangeError("UserError: unexpected '" + String.fromCharCode(c) + "', '" + String.fromCharCode(characterCode) + "' expected");
    }
    return Input.trimLeft(input, position, 1);
  };
  Input.getFirst = function (input, position) {
    return position < input.length ? input.charCodeAt(position) : Input.EOF;
  };
  Input.startsWith = function (input, position, s) {
    var length = s.length;
    if (position + length > input.length) {
      return false;
    }
    var i = -1;
    while (++i < length) {
      if (asciiToLowerCase(input.charCodeAt(position + i)) !== asciiToLowerCase(s.charCodeAt(i))) {
        return false;
      }
    }
    return true;
  };
  Input.exec = function (input, position, regularExpression) {
    if (position === input.length) {
      return undefined;
    }
    var match = regularExpression.exec(position === 0 ? input : input.slice(position));
    return match != undefined ? match[0] : undefined;
  };

  function ParseResult(result, position) {
    this.result = result;
    this.position = position;
  }

  var parseMatrix = function (input, position, context) {
    var openingBracket = "{".charCodeAt(0);
    var closingBracket = "}".charCodeAt(0);

    //position = Input.parseCharacter(input, position, openingBracket);
    var rows = [];
    var firstRow = true;
    while (firstRow || Input.getFirst(input, position) === ",".charCodeAt(0)) {
      if (firstRow) {
        firstRow = false;
      } else {
        position = Input.trimLeft(input, position, 1);
      }
      position = Input.parseCharacter(input, position, openingBracket);
      var row = [];
      var firstCell = true;
      while (firstCell || Input.getFirst(input, position) === ",".charCodeAt(0)) {
        if (firstCell) {
          firstCell = false;
        } else {
          position = Input.trimLeft(input, position, 1);
        }
        var tmp = parseExpression(input, position, context, true, 0, undefined);
        position = tmp.position;
        row.push(tmp.result);
      }
      position = Input.parseCharacter(input, position, closingBracket);
      rows.push(row);
    }
    //position = Input.parseCharacter(input, position, closingBracket);
    return new ParseResult(context.wrap(Expression.Matrix.fromArray(rows)), position);
  };

  var parseLaTeXMatrix = function (input, position, context, kind) {
    position = Input.trimLeft(input, position, ("\\begin{" + kind + "}").length);
    var rows = [];
    var firstRow = true;
    while (firstRow || (Input.getFirst(input, position) === "\\".charCodeAt(0) && Input.startsWith(input, position, "\\\\"))) {
      if (firstRow) {
        firstRow = false;
      } else {
        position = Input.trimLeft(input, position, 2);
      }
      var row = [];
      var firstCell = true;
      while (firstCell || Input.getFirst(input, position) === "&".charCodeAt(0)) {
        if (firstCell) {
          firstCell = false;
        } else {
          position = Input.trimLeft(input, position, 1);
        }
        var tmp = parseExpression(input, position, context, true, ADDITIVE_PRECEDENCE - 1, undefined);
        position = tmp.position;
        row.push(tmp.result);
      }
      rows.push(row);
    }
    if (Input.startsWith(input, position, "\\end{" + kind + "}")) {
      position = Input.trimLeft(input, position, ("\\end{" + kind + "}").length);
    }
    return new ParseResult(context.wrap(Expression.Matrix.fromArray(rows)), position);
  };

  var parseLaTeXArgument = function (input, position, context) {
    return parseExpression(input, position, context, true, 0, undefined);
  };

  var getDecimalFraction = function (integerPartAsString, nonRepeatingFractionalPartAsString, repeatingFractionalPartAsString, exponentSingPartAsString, exponentPartAsString) {
    var numerator = Expression.ZERO;
    var denominator = undefined;
    var factor = undefined;

    if (integerPartAsString != undefined) {
      numerator = Expression.Integer.fromString(integerPartAsString);
    }
    if (nonRepeatingFractionalPartAsString != undefined) {
      factor = Expression.pow(Expression.TEN, nonRepeatingFractionalPartAsString.length);
      numerator = numerator.multiply(factor).add(Expression.Integer.fromString(nonRepeatingFractionalPartAsString));
      denominator = denominator == undefined ? factor : denominator.multiply(factor);
    }
    if (repeatingFractionalPartAsString != undefined) {
      factor = Expression.pow(Expression.TEN, repeatingFractionalPartAsString.length).subtract(Expression.ONE);
      numerator = numerator.multiply(factor).add(Expression.Integer.fromString(repeatingFractionalPartAsString));
      denominator = denominator == undefined ? factor : denominator.multiply(factor);
    }
    if (exponentPartAsString != undefined) {
      factor = Expression.pow(Expression.TEN, Number.parseInt(exponentPartAsString, 10));
      if (exponentSingPartAsString === "-") {
        denominator = denominator == undefined ? factor : denominator.multiply(factor);
      } else {
        numerator = numerator.multiply(factor);
      }
    }

    var value = denominator == undefined ? numerator : numerator.divide(denominator);
    return value;
  };

  var parseDecimalFraction = function (input, position, context, isMatrixElement) {
    var isOnlyInteger = true;
    var match = undefined;
    var integerPartAsString = undefined;
    var nonRepeatingFractionalPartAsString = undefined;
    var repeatingFractionalPartAsString = undefined;
    if ((match = Input.exec(input, position, digits)) != undefined) {
      integerPartAsString = match;
      position += match.length;
    }
    var c = Input.getFirst(input, position);
    if (c === ".".charCodeAt(0) || (!isMatrixElement && c === ",".charCodeAt(0))) {
      isOnlyInteger = false;
      position += 1;
      if ((match = Input.exec(input, position, digits)) != undefined) {
        nonRepeatingFractionalPartAsString = match;
        position += match.length;
      }
      c = Input.getFirst(input, position);
      if (c === "(".charCodeAt(0)) {
        if ((match = Input.exec(input, position + 1, digits)) != undefined) {
          c = Input.getFirst(input, position + 1 + match.length);
          if (c === ")".charCodeAt(0)) {
            position += 1 + match.length + 1;
            repeatingFractionalPartAsString = match;
          }
        }
      }
    }
    var result = undefined;
    if (integerPartAsString != undefined || nonRepeatingFractionalPartAsString != undefined || repeatingFractionalPartAsString != undefined) {
      var exponentSingPartAsString = "";
      var exponentPartAsString = undefined;
      c = Input.getFirst(input, position);
      if (c === "e".charCodeAt(0) || c === "E".charCodeAt(0)) {
        c = Input.getFirst(input, position + 1);
        if (c === "+".charCodeAt(0) || c === "-".charCodeAt(0)) {
          exponentSingPartAsString = input.slice(position + 1, position + 2);
        }
        if ((match = Input.exec(input, position + 1 + exponentSingPartAsString.length, digits)) != undefined) {
          position += 1 + exponentSingPartAsString.length + match.length;
          exponentPartAsString = match;
        }
      }
      result = getDecimalFraction(integerPartAsString, nonRepeatingFractionalPartAsString, repeatingFractionalPartAsString, exponentSingPartAsString, exponentPartAsString);
      result = context.wrap(result);
      position = Input.trimLeft(input, position, 0);
    }
    //!
    if (isOnlyInteger || result == undefined) {
      if ((match = Input.exec(input, position, vulgarFractions)) != undefined) {
        var tmp = parseExpression(normalizeVulgarFractions(match), 0, context, isMatrixElement, 0, undefined);
        if (result != undefined) {
          result = ADDITION.i(result, tmp.result).addPosition(position, ADDITION.name.length, input);
        } else {
          result = tmp.result;
        }
        position = Input.trimLeft(input, position, match.length);
      }
    }
    return result != undefined ? new ParseResult(result, position) : undefined;
  };

  // TODO: sticky flags - /\s+/y
  var whiteSpaces = /^\s+/;
  var digits = /^\d+/;
  // Base Latin, Base Latin upper case, Base Cyrillic, Base Cyrillic upper case, Greek alphabet
  var symbols = /^(?:alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|mu|nu|xi|omicron|pi|rho|varsigma|sigma|tau|upsilon|phi|chi|psi|omega|[a-zA-Z\u0430-\u044F\u0410-\u042F\u03B1-\u03C9])(?:\_\d+|\_\([a-z\d]+,[a-z\d]+\)|[\u2080-\u2089]+)?/;
  var superscripts = /^[\u00B2\u00B3\u00B9\u2070\u2074-\u2079]+/;
  var vulgarFractions = /^[\u00BC-\u00BE\u2150-\u215E]/;


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
      if (charCode === 0x2070) {
        return "0";
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

  var normalizeFullwidthForms = function (charCode) {
    return String.fromCharCode(charCode - 0xFF01 + 0x0021);
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

  var parseExpression = function (input, position, context, isMatrixElement, precedence, left) {
    var ok = true;
    var firstCharacterCode = Input.getFirst(input, position);
    var isDecimalFraction = false;
    //!

    while (firstCharacterCode !== Input.EOF && ok) {
      var op = undefined;
      var operand = undefined;
      var tmp = undefined;
      var match = undefined;
      var right = undefined;

      if (firstCharacterCode < "0".charCodeAt(0) || firstCharacterCode > "9".charCodeAt(0)) {
        var operationsArray = operationSearchCache.getAll(firstCharacterCode);
        if (operationsArray != undefined) {
          var length = operationsArray.length;
          var bestMatchLength = 0;//? "^T" and "^"
          var j = -1;
          while (++j < length) {
            var candidate = operationsArray[j];
            if ((left != undefined && (candidate.arity !== 1 || candidate.rightToLeftAssociative !== RIGHT_TO_LEFT || precedence < MULTIPLICATION.precedence) ||
                 left == undefined && candidate.arity === 1 && candidate.rightToLeftAssociative === RIGHT_TO_LEFT) &&
                Input.startsWith(input, position, candidate.name) &&
                //(!candidate.xyz || !isAlpha(Input.getFirst(input, position + candidate.name.length))) &&//TODO: fix - ExpressionParser.parse("George")
                bestMatchLength < candidate.name.length) {
              op = candidate;
              bestMatchLength = op.name.length;
            }
          }
        }
      }
      if (op != null && op.name === "\\" && isMatrixElement) {//TODO: optimize
        op = null;
      }

      if (firstCharacterCode === "f".charCodeAt(0) && op == undefined && (left == undefined && precedence <= UNARY_PRECEDENCE || precedence < MULTIPLICATION.precedence) && Input.startsWith(input, position, "frac")) { // !isAlpha(Input.getFirst(input, position + "frac".length))
        // https://en.wikipedia.org/wiki/Operand#Positioning_of_operands - prefix notation
        position = Input.trimLeft(input, position, "frac".length);
        tmp = parseExpression(input, position, context, isMatrixElement, MULTIPLICATION.precedence, undefined);
        var a = tmp.result;
        position = tmp.position;
        tmp = parseExpression(input, position, context, isMatrixElement, MULTIPLICATION.precedence, undefined);
        var b = tmp.result;
        position = tmp.position;
        // addPosition - ?
        operand = a.divide(b);
        ok = true;
      } else if (op != undefined) {
        if (precedence > op.precedence + (op.rightToLeftAssociative === RIGHT_TO_LEFT ? 0 : -1)) {
          ok = false;
        } else {
          var operatorPosition = position;
          position = Input.trimLeft(input, position, op.name.length);
          if (op.arity === 1 && op.rightToLeftAssociative !== RIGHT_TO_LEFT) {
            //TODO: fix
            ExpressionParser.startPosition = operatorPosition;
            ExpressionParser.endPosition = operatorPosition + op.name.length;
            ExpressionParser.input = input;
            left = op.i(left).addPosition(operatorPosition, op.name.length, input);
          } else {
            if (op.arity === 1 && op.rightToLeftAssociative === RIGHT_TO_LEFT && op.precedence === UNARY_PRECEDENCE_PLUS_ONE && op.name.length > 1 &&
                (op.name === "sin" || op.name === "cos" || op.name === "sen" || op.name === "tan" || op.name === "tg") &&
                (Input.startsWith(input, position, EXPONENTIATION.name) || (match = Input.exec(input, position, superscripts)) != null)) {
              
              // cos^2(x)
              //!new 2017-11-04
              // parse an operator for the exponentiation
              var exponentiationPosition = position;

              var exponentiationLength = 0;
              var middle = null;
              if (match != null) {
                exponentiationLength = match.length;
                position = position + exponentiationLength;
                middle = Expression.Integer.fromString(normalizeSuperscripts(match));
              } else {
                exponentiationLength = EXPONENTIATION.name.length;
                position = position + exponentiationLength;
                tmp = parseExpression(input, position, context, isMatrixElement, EXPONENTIATION.precedence, undefined);
                middle = tmp.result;
                position = tmp.position;
              }

              // parse an operator for the current operator
              tmp = parseExpression(input, position, context, isMatrixElement, op.precedence, undefined);
              right = tmp.result;
              position = tmp.position;
              operand = EXPONENTIATION.i(op.i(right).addPosition(operatorPosition, op.name.length, input), middle).addPosition(exponentiationPosition, exponentiationLength, input);
            } else {
              tmp = parseExpression(input, position, context, isMatrixElement, op.precedence, undefined);
              right = tmp.result;
              position = tmp.position;
              //TODO: fix `1/(2-2)`
              ExpressionParser.startPosition = operatorPosition;
              ExpressionParser.endPosition = operatorPosition + op.name.length;
              ExpressionParser.input = input;
              if (op.arity === 1) {
                // left <implicit multiplication> operand
                operand = op.i(right).addPosition(operatorPosition, op.name.length, input);
              } else if (op.arity === 2) {
                left = op.i(left, right).addPosition(operatorPosition, op.name.length, input);
              } else {
                throw new RangeError();
              }
            }
          }
        }
      } else if (left == undefined || precedence < MULTIPLICATION.precedence || (precedence === UNARY_PRECEDENCE_PLUS_ONE && isDecimalFraction && (match = Input.exec(input, position, symbols)) != undefined)) {
        if (firstCharacterCode === "(".charCodeAt(0)) {
          position = Input.parseCharacter(input, position, "(".charCodeAt(0));
          tmp = parseExpression(input, position, context, false, 0, undefined);
          operand = tmp.result;
          position = tmp.position;
          position = Input.parseCharacter(input, position, ")".charCodeAt(0));
        } else if (firstCharacterCode === "{".charCodeAt(0)) {
          position = Input.parseCharacter(input, position, "{".charCodeAt(0));
          if (Input.getFirst(input, position) === "{".charCodeAt(0)) {
            tmp = parseMatrix(input, position, context);
            operand = tmp.result;
            position = tmp.position;
          } else {
            tmp = parseLaTeXArgument(input, position, context);
            operand = tmp.result;
            position = tmp.position;
          }
          position = Input.parseCharacter(input, position, "}".charCodeAt(0));
        } else if ((tmp = parseDecimalFraction(input, position, context, isMatrixElement)) != undefined) {
          operand = tmp.result;
          position = tmp.position;
          isDecimalFraction = true;
        } else if (firstCharacterCode === "\\".charCodeAt(0) && (Input.startsWith(input, position, "\\begin{bmatrix}") || Input.startsWith(input, position, "\\begin{vmatrix}") || Input.startsWith(input, position, "\\begin{pmatrix}") || Input.startsWith(input, position, "\\begin{matrix}"))) {
          var kind = "matrix";
          if (Input.startsWith(input, position, "\\begin{bmatrix}")) {
            kind = "bmatrix";
          } else if (Input.startsWith(input, position, "\\begin{vmatrix}")) {
            kind = "vmatrix";
          } else if (Input.startsWith(input, position, "\\begin{pmatrix}")) {
            kind = "pmatrix";
          }
          tmp = parseLaTeXMatrix(input, position, context, kind);
          operand = tmp.result;
          if (Input.startsWith(input, position, "\\begin{vmatrix}")) {
            operand = operand.determinant();//!
          }
          position = tmp.position;
        } else if ((match = Input.exec(input, position, symbols)) != undefined) {
          var symbolName = match;
          symbolName = normalizeSubscripts(symbolName);
          symbolName = normalizeGreek(symbolName);
          operand = context.get(symbolName);
          if (operand == undefined) {
            //TODO: move to context - ?
            operand = symbolName === "I" || symbolName === "U" || symbolName === "E" ? new Expression.IdentityMatrix(symbolName) : new Expression.Symbol(symbolName);
            operand = context.wrap(operand);
          } else {
            operand = context.wrap(operand);
          }
          position = Input.trimLeft(input, position, match.length);
        } else if (firstCharacterCode === "|".charCodeAt(0) && left == undefined) {
          position = Input.parseCharacter(input, position, "|".charCodeAt(0));
          tmp = parseExpression(input, position, context, isMatrixElement, 0, undefined);
          operand = tmp.result;
          position = tmp.position;
          position = Input.parseCharacter(input, position, "|".charCodeAt(0));
          
          operand = operand.determinant();//!
        } else {
          ok = false;
        }
      } else {
        ok = false;
      }

      //!TODO: fix
      if (!ok && left != undefined && precedence <= EXPONENTIATION.precedence + (EXPONENTIATION.rightToLeftAssociative === RIGHT_TO_LEFT ? 0 : -1)) {
        if ((match = Input.exec(input, position, superscripts)) != undefined) {
          // implicit exponentiation
          //TODO: check position
          left = EXPONENTIATION.i(left, Expression.Integer.fromString(normalizeSuperscripts(match))).addPosition(position, EXPONENTIATION.name.length, input);
          position = Input.trimLeft(input, position, match.length);
          ok = true;//!
        }
      }

      if (!ok && firstCharacterCode === "\\".charCodeAt(0) && !Input.startsWith(input, position, "\\\\")) { // isAlpha(Input.getFirst(input, position + 1))
        if (!Input.startsWith(input, position + 1, "end{bmatrix}") && !Input.startsWith(input, position + 1, "end{vmatrix}") && !Input.startsWith(input, position + 1, "end{pmatrix}") && !Input.startsWith(input, position + 1, "end{matrix}")) {
        // TODO: LaTeX - ?
        ok = true;
        position += 1;
        }
      }

      if (operand != undefined) {
        if (left != undefined) {
          // implied multiplication
          tmp = parseExpression(input, position, context, isMatrixElement, MULTIPLICATION.precedence, operand);
          var right1 = tmp.result;
          var position1 = tmp.position;
          left = MULTIPLICATION.i(left, right1).addPosition(position, MULTIPLICATION.name.length, input);
          position = position1;
        } else {
          left = operand;
        }
      }
      firstCharacterCode = Input.getFirst(input, position);
    }

    if (left == undefined) {
      ExpressionParser.startPosition = position;
      ExpressionParser.endPosition = position + 1;
      ExpressionParser.input = input;
      var c = Input.getFirst(input, position);
      if (c === Input.EOF) {
        throw new RangeError("UserError: unexpected end of input");//TODO: fix
      }
      throw new RangeError("UserError: unexpected '" + String.fromCharCode(c) + "'");//TODO: fix
    }
    return new ParseResult(left, position);
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
      return {code: code - 0x0660 + 0x0030, name: "arab"};
    }
    if (code >= 0x06F0 && code <= 0x06F9) {
      return {code: code - 0x06F0 + 0x0030, name: "arabext"};
    }
    if (code >= 0x0966 && code <= 0x096F) {
      return {code: code - 0x0966 + 0x0030, name: "deva"};
    }
    if (code >= 0x09E6 && code <= 0x09EF) {
      return {code: code - 0x09E6 + 0x0030, name: "beng"};
    }
    if (code >= 0x0A66 && code <= 0x0A6F) {
      return {code: code - 0x0A66 + 0x0030, name: "guru"};
    }
    if (code >= 0x0AE6 && code <= 0x0AEF) {
      return {code: code - 0x0AE6 + 0x0030, name: "gujr"};
    }
    if (code >= 0x0B66 && code <= 0x0B6F) {
      return {code: code - 0x0B66 + 0x0030, name: "orya"};
    }
    if (code >= 0x0BE6 && code <= 0x0BEF) {
      return {code: code - 0x0BE6 + 0x0030, name: "tamldec"};
    }
    if (code >= 0x0C66 && code <= 0x0C6F) {
      return {code: code - 0x0C66 + 0x0030, name: "telu"};
    }
    if (code >= 0x0CE6 && code <= 0x0CEF) {
      return {code: code - 0x0CE6 + 0x0030, name: "knda"};
    }
    if (code >= 0x0D66 && code <= 0x0D6F) {
      return {code: code - 0x0D66 + 0x0030, name: "mlym"};
    }
    if (code >= 0x0E50 && code <= 0x0E59) {
      return {code: code - 0x0E50 + 0x0030, name: "thai"};
    }
    if (code >= 0x0ED0 && code <= 0x0ED9) {
      return {code: code - 0x0ED0 + 0x0030, name: "laoo"};
    }
    if (code >= 0x0F20 && code <= 0x0F29) {
      return {code: code - 0x0F20 + 0x0030, name: "tibt"};
    }
    if (code >= 0x1040 && code <= 0x1049) {
      return {code: code - 0x1040 + 0x0030, name: "mymr"};
    }
    if (code >= 0x17E0 && code <= 0x17E9) {
      return {code: code - 0x17E0 + 0x0030, name: "khmr"};
    }
    if (code >= 0x1810 && code <= 0x1819) {
      return {code: code - 0x1810 + 0x0030, name: "mong"};
    }
    if (code >= 0x1946 && code <= 0x194F) {
      return {code: code - 0x1946 + 0x0030, name: "limb"};
    }
    if (code >= 0x1B50 && code <= 0x1B59) {
      return {code: code - 0x1B50 + 0x0030, name: "bali"};
    }
    if (code >= 0xFF10 && code <= 0xFF19) {
      return {code: code - 0xFF10 + 0x0030, name: "fullwide"};
    }
    var c = replaceHanidec(code);
    if (c !== -1) {
      return {code: c, name: "hanidec"};
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
      return normalizeFullwidthForms(charCode);
    }
    if (charCode === 0x060C || charCode === 0x066B) {
      return ",";
    }
    var y = replaceSimpleDigit(charCode);
    if (y != undefined) {
      // TODO: remove
      if (typeof hit === "function") {
        hit({digit: y.name});
      }
      return String.fromCharCode(y.code);
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
    return undefined;
  };
  //var replaceRegExp = /[...]/g;
  //var replaceFunction = function (c) {
  //  return getCharCodeReplacement(c.charCodeAt(0));
  //};
  //input = input.replace(replaceRegExp, replaceFunction); - slow in Chrome
  var replaceSomeChars = function (input) {
    var lastIndex = 0;
    var result = "";
    for (var i = 0; i < input.length; i += 1) {
      var charCode = input.charCodeAt(i);
      if (charCode > 0x007F || charCode === 0x003A) {
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

    var position = 0;
    position = Input.trimLeft(input, position, 0);
    var tmp = parseExpression(input, position, context, false, 0, undefined);
    var c = Input.getFirst(input, tmp.position);
    if (c !== Input.EOF) {
      ExpressionParser.startPosition = tmp.position;
      ExpressionParser.endPosition = tmp.position + 1;
      ExpressionParser.input = input;
      throw new RangeError("UserError: unexpected '" + String.fromCharCode(c) + "'");
    }

    return tmp.result;
  };

  ExpressionParser.startPosition = -1;
  ExpressionParser.endPosition = -1;
  ExpressionParser.input = "";

  var getConstant = function (e) {
    if (e === "pi" || e === "\u03C0") {
      return Expression.PI;
    }
    if (e === "e") {
      return Expression.E;
    }
    if (e === "i") {
      return Expression.I;
    }
    return undefined;
  };

  ExpressionParser.Context = function (getter, needsWrap) {
    this.getter = getter;
    this.needsWrap = needsWrap == undefined ? true : needsWrap;
  };
  ExpressionParser.Context.prototype.get = function (e) {
    if (this.getter != undefined) {
      var x = this.getter(e);
      if (x != undefined) {
        return x;
      }
    }
    return getConstant(e);
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
    operationSearchCache.append(newOperation.name.charCodeAt(0), newOperation);
  };

  ExpressionParser.addDenotations = function (operationName, denotations) {
    var os = operationSearchCache.getAll(operationName.charCodeAt(0));
    var operation = undefined;
    var i = -1;
    while (++i < os.length) {
      var o = os[i];
      if (o.name === operationName) {
        operation = o;
      }
    }
    var added = {};
    added[operationName] = true;
    for (var key in denotations) {
      if (Object.prototype.hasOwnProperty.call(denotations, key)) {
        var denotation = denotations[key];
        if (added[denotation] == undefined) {
          added[denotation] = true;
          var newOperation = new Operator(denotation, operation.arity, operation.rightToLeftAssociative, operation.precedence, operation.i);
          //operations.push(newOperation);
          operationSearchCache.append(newOperation.name.charCodeAt(0), newOperation);
        }
      }
    }
  };

  export default ExpressionParser;
