
const baseURL = new URL('file://');
baseURL.pathname = `${process.cwd()}/`;

export async function resolve(specifier,
                              parentModuleURL = baseURL,
                              defaultResolver) {
  if (!specifier.endsWith('.js')) {
    return defaultResolver(specifier, parentModuleURL);
  }
  return {
    url: new URL(specifier, parentModuleURL).href,
    format: 'esm'
  };
}
