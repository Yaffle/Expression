/*jshint esversion:6, bitwise:false*/

// https://webassembly.github.io/wabt/demo/wat2wasm/index.html helps a lot

function wast2wasm(sExpression) {

  function sExpressionToJSON(s) {
    return JSON.parse(s.trim().replace(/\s*\)/g, ']').replace(/\(\s*/g, '[').replace(/\s+/g, ',').replaceAll('"', '').replace(/([^\[\],]+)/g, '"$1"'));
  }

  const bytes = [];

  function pushByte(byte) {
    if ((byte & 0x7F) !== byte) {
      throw new RangeError();
    }
    bytes.push(byte);
  }

  function pushUnsignedInteger(position, n) {
    do {
      const byte = (n & 0x7F) | ((n >>> 7) !== 0 ? 0x80 : 0x00);
      n >>>= 7;
      bytes.splice(position, 0, byte);
      position += 1;
    } while (n !== 0);
  }

  function pushSize(f) {
    const start = bytes.length;
    f();
    const end = bytes.length;
    pushUnsignedInteger(start, end - start);
  }

  function pushString(s) {
    pushByte(s.length);
    for (let i = 0; i < s.length; i++) {
      pushByte(s.charCodeAt(i));
    }
  }

  function pushOpcode(opcode) {
    if (opcode <= 0xF0) {
      bytes.push(opcode);
    } else {
      bytes.push((opcode >> 8));
      pushUnsignedInteger(bytes.length, opcode & 0xFF);
    }
  }

  function magic() {
    for (const byte of [0x00, 0x61, 0x73, 0x6d]) {
      pushByte(byte);
    }
  }

  function version() {
    for (const byte of [0x01, 0x00, 0x00, 0x00]) {
      pushByte(byte);
    }
  }

  function emitTypes(node, expectedName) {
    console.assert(node[0] === expectedName);
    pushSize(function () {
      for (let j = 1; j < node.length; j++) {
        pushByte(type(node[j]));
      }
    });
  }

  function typesec() {
    pushByte(0x01);
    pushSize(function () {
      pushByte(types.length);
      for (const typeNode of types) {
        const node = typeNode[2];
        if (node[0] === 'func') {
          pushByte(0x60);
          emitTypes(node[1], 'param');
          emitTypes(node[2], 'result');
        } else {
          throw new TypeError();
        }
      }
    });
  }

  function importsec() {
    pushByte(0x02);
    pushSize(function () {
      pushByte(imports.length);
      for (const importNode of imports) {
        // (import "env" "memory" (memory $0 0))
        console.assert(importNode.length === 4);
        pushString(importNode[1]);
        pushString(importNode[2]);
        if (importNode[3][0] === 'memory') {
          pushByte(0x02);
          pushByte(0x00);
          pushByte(0x00);
        } else {
          throw new TypeError();
        }
      }
    });
  }

  function funcsec() {
    pushByte(0x03);
    pushSize(function () {
      pushByte(funcs.length);
      for (let i = 0; i < funcs.length; i++) {
        pushByte(i);
      }
    });
  }

  function exportsec() {
    pushByte(0x07);
    pushSize(function () {
      pushByte(exports.length);
      // (export "add" (func $module/add))
      for (const exportNode of exports) {
        console.assert(exportNode.length === 3);
        pushString(exportNode[1]);
        if (exportNode[2][0] === 'func') {
          console.assert(exportNode[2].length === 2);
          pushByte(0x00);
          const internalName = exportNode[2][1];
          const index = funcs.findIndex(func => func[1] === internalName);
          pushByte(index);
        } else {
          throw new TypeError();
        }
      }
    });
  }

  const Types = {
    'i32': 0x7F,
    'i64': 0x7E,
    'f32': 0x7D,
    'f64': 0x7C,
    'v128': 0x7B
  };

  function type(name) {
    const code = Types[name];
    if (code == undefined) {
      throw new TypeError('unknown type: ' + name);
    }
    return code;
  }

  // from https://pengowray.github.io/wasm-ops/:
  const OpcodesList = 'unreachable nop block loop if else - - - - - end br br_if br_table return call call_indirect - - - - - - - - drop select - - - - local.get local.set local.tee global.get global.set - - - i32.load i64.load f32.load f64.load i32.load8_s i32.load8_u i32.load16_s i32.load16_u i64.load8_s i64.load8_u i64.load16_s i64.load16_u i64.load32_s i64.load32_u i32.store i64.store f32.store f64.store i32.store8 i32.store16 i64.store8 i64.store16 i64.store32 memory.size memory.grow i32.const i64.const f32.const f64.const i32.eqz i32.eq i32.ne i32.lt_s i32.lt_u i32.gt_s i32.gt_u i32.le_s i32.le_u i32.ge_s i32.ge_u i64.eqz i64.eq i64.ne i64.lt_s i64.lt_u i64.gt_s i64.gt_u i64.le_s i64.le_u i64.ge_s i64.ge_u f32.eq f32.ne f32.lt f32.gt f32.le f32.ge f64.eq f64.ne f64.lt f64.gt f64.le f64.ge i32.clz i32.ctz i32.popcnt i32.add i32.sub i32.mul i32.div_s i32.div_u i32.rem_s i32.rem_u i32.and i32.or i32.xor i32.shl i32.shr_s i32.shr_u i32.rotl i32.rotr i64.clz i64.ctz i64.popcnt i64.add i64.sub i64.mul i64.div_s i64.div_u i64.rem_s i64.rem_u i64.and i64.or i64.xor i64.shl i64.shr_s i64.shr_u i64.rotl i64.rotr f32.abs f32.neg f32.ceil f32.floor f32.trunc f32.nearest f32.sqrt f32.add f32.sub f32.mul f32.div f32.min f32.max f32.copysign f64.abs f64.neg f64.ceil f64.floor f64.trunc f64.nearest f64.sqrt f64.add f64.sub f64.mul f64.div f64.min f64.max f64.copysign i32.wrap_i64 i32.trunc_f32_s i32.trunc_f32_u i32.trunc_f64_s i32.trunc_f64_u i64.extend_i32_s i64.extend_i32_u i64.trunc_f32_s i64.trunc_f32_u i64.trunc_f64_s i64.trunc_f64_u f32.convert_i32_s f32.convert_i32_u f32.convert_i64_s f32.convert_i64_u f32.demote_f64 f64.convert_i32_s f64.convert_i32_u f64.convert_i64_s f64.convert_i64_u f64.promote_f32 i32.reinterpret_f32 i64.reinterpret_f64 f32.reinterpret_i32 f64.reinterpret_i64'.split(' ');

  const Opcodes = {};
  for (let i = 0; i < OpcodesList.length; i++) {
    Opcodes[OpcodesList[i]] = i;
  }

  const SIMDOpcodes = {
    'i8x16.splat': 0xFD0F,
    'i8x16.ge_u': 0xFD2C,
    'i8x16.all_true': 0xFD63,

    'i32x4.splat': 0xFD11,
    'i32x4.lt_s': 0xFD39,
    'i32x4.ge_s': 0xFD3F,
    'i32x4.gt_s': 0xFD3B,
    'i32x4.all_true': 0xFDA3,
    'v128.any_true': 0xFD53,
    'i8x16.gt_s': 0xFD27,
    'i8x16.gt_u': 0xFD28,

    'v128.xor': 0xFD51,
    'v128.load': 0xFD00,
    'v128.store': 0xFD0B
  };

  const AllOpcodes = Object.assign({}, Opcodes, SIMDOpcodes);

  function opcode(name) {
    const code = AllOpcodes[name];
    if (code == undefined) {
      throw new TypeError('unknown instruction: ' + name);
    }
    return code;
  }

  const voidBlockType = 0x40;

  const labels = [];
  let vars = [];

  function localInstr(node) {
    if (node.length > 3) {
      throw new TypeError();
    }
    const variableName = node[1];
    const index = vars.indexOf(variableName);
    if (node.length === 3) {
      emitCode(node[2]);
    }
    pushByte(opcode(node[0]));
    pushByte(index);
  }

  function constInstr(node) {
    if (node.length !== 2) {
      throw new TypeError();
    }
    //! 0 <= c <= 63
    if ((Number(node[1]) & 0x3F).toString() !== node[1]) {
      throw new RangeError('unsupported literal: ' + node[1]);
    }
    pushByte(opcode(node[0]));
    pushByte(Number(node[1]));
  }

  function blockInstr(node, withIf) {
    const label = typeof node[withIf ? 2 : 1] === 'string' ? node[withIf ? 2 : 1] : null;
    labels.push(label);
    if (withIf) {
      emitCode(node[1]);
    }
    pushByte(opcode(node[0]));
    pushByte(voidBlockType);
    const s = 1 + (label != null ? 1 : 0) + (withIf ? 1 : 0);
    if (withIf) {
      if (node.length > s + 2) {
        throw new TypeError();
      }
      emitCode(node[s]);
      if (node.length === s + 2) {
        pushByte(opcode('else'));
        emitCode(node[s + 1]);
      }
    } else {
      for (let i = s; i < node.length; i++) {
        emitCode(node[i]);
      }
    }
    pushByte(opcode('end'));
    labels.pop();
  }

  function brInstr(node, withIf) {
    if (node.length !== 2 + (withIf ? 1 : 0)) {
      throw new TypeError();
    }
    const label = node[1];
    if (withIf) {
      emitCode(node[2]);
    }
    pushByte(opcode(node[0]));
    pushByte(labels.length - 1 - labels.lastIndexOf(label));
  }

  function memoryInstr(node) {
    const tmp1 = /^[ifv](32|64|128)\.(?:store|load)(8|16|32|64|)(?:_u|_s)?$/.exec(node[0])
    if (tmp1 == null) {
      throw new TypeError(node[0]);
    }
    const alignment = Math.round(Math.log2(Number(tmp1[2] || tmp1[1]) / 8));
    let offset = 0;
    for (let i = 1; i < node.length; i++) {
      if (typeof node[i] === 'string') {
        const tmp = /^offset=(\d+)$/.exec(node[i]);
        if (tmp == null) {
          throw new TypeError();
        }
        offset = Number(tmp[1]);
      } else {
        emitCode(node[i]);
      }
    }
    pushOpcode(opcode(node[0]));
    pushByte(alignment);
    pushByte(offset);
  }

  function callInstr(node) {
    if (node.length < 2) {
      throw new TypeError();
    }
    const functionName = node[1];
    const index = funcs.findIndex(func => func[1] === functionName);
    for (let i = 2; i < node.length; i++) {
      emitCode(node[i]);
    }
    pushByte(opcode(node[0]));
    pushByte(index);
  }

  function emitCode(node) {
    if (node[0].startsWith('local.')) {
      localInstr(node);
    } else if (node[0].endsWith('.const')) {
      constInstr(node);
    } else if (node[0] === 'loop' || node[0] === 'block' || node[0] === 'if') {
      blockInstr(node, node[0] === 'if');
    } else if (node[0] === 'br' || node[0] === 'br_if') {
      brInstr(node, node[0] === 'br_if');
    } else if (node[0].indexOf('.store') !== -1 || node[0].indexOf('.load') !== -1) {
      memoryInstr(node);
    } else if (node[0] === 'call') {
      callInstr(node);
    } else {
      for (let i = 1; i < node.length; i++) {
        emitCode(node[i]);
      }
      pushOpcode(opcode(node[0]));
    }
  }

  function codesec() {
    pushByte(0x0A);
    pushSize(function () {
      pushByte(funcs.length);
      // (func $module/add (param $a i32) (param $b i32) (result i32)
      // (return (i32.add (local.get $a) (local.get $b))))
      for (const func of funcs) {
        pushSize(function () {
          // skip func[1]
          // skip param and result
          const locals = [];
          const params = [];
          const instrs = [];
          for (let i = 2; i < func.length; i++) {
            if (func[i][0] === 'param') {
              params.push({name: func[i][1], type: func[i][2]});
            } else if (func[i][0] === 'result') {
              
            } else if (func[i][0] === 'local') {
              console.assert(func[i].length === 3);
              locals.push({name: func[i][1], type: func[i][2]});
            } else {
              instrs.push(func[i]);
            }
          }
          pushByte(locals.length);
          for (const local of locals) {
            pushByte(0x01);
            pushByte(type(local.type));
          }
          vars = params.concat(locals).map(x => x.name);
          for (const instr of instrs) {
            emitCode(instr);
          }
          pushByte(opcode('end'));
        });
      }
    });
  }

  const types = [];
  const imports = [];
  const funcs = [];
  const exports = [];

  const moduleNode = sExpressionToJSON(sExpression);
  if (moduleNode[0] !== 'module') {
    throw new TypeError();
  }
  for (let i = 1; i < moduleNode.length; i++) {
    const node = moduleNode[i];
    if (node[0] === 'type') {
      types.push(node);
    } else if (node[0] === 'import') {
      imports.push(node);
    } else if (node[0] === 'func') {
      funcs.push(node);
    } else if (node[0] === 'export') {
      exports.push(node);
    } else {
      throw new TypeError();
    }
  }

  magic();
  version();
  typesec();
  importsec();
  funcsec();
  //tablesec();
  //memsec();
  //globalsec();
  exportsec();
  //startsec();
  //elemsec();
  //datacountsec();
  codesec();
  //datasec();
  
  return new Uint8Array(bytes);
}

export default wast2wasm;