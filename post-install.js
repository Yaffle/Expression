// post-install.js

/**
 * Script to run after npm install
 *
 * Make a symbolic lync to node_modules
 */


import fs from 'fs';

fs.symlink('../../../node_modules', './node_modules', 'junction', (err) => { console.log(err); });
