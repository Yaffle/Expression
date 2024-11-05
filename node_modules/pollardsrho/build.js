
const fs = require('fs');

const name = process.argv[2] || '';

const filePath = './' + name + '.wasm';
const filePath1 = './' + name + '.js';

fs.readFile(filePath, (err, data) => {
  if (err) {
    console.error('Error reading file:', err);
    return;
  }
  const byteArr = Array.from(data);
  const textString = `
    let ${name} = typeof WebAssembly !== 'undefined' ? new WebAssembly.Instance(new WebAssembly.Module(new Uint8Array([${byteArr.join(', ')}]))).exports.${name} : null;
    export default ${name};
  `;
  fs.writeFile(filePath1, textString, (err) => {
    if (err) {
      console.error('Error writing to file:', err);
      return;
    }
  });
});
