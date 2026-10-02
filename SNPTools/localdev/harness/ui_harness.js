/* ui_harness.js -- load SNPTools index.html in jsdom exactly as a browser would
 * (all <script>s, same order), with:
 *   - static files served from <siteRoot>
 *   - *.php requests executed by the real PHP scripts through the php CLI
 *     (php_shim.php), so processForm.php -> h5_to_vcf.py -> HDF5 runs for real
 *   - d3 v5 served from node_modules (the page loads it from d3js.org)
 * Exposes openSite() -> {window, $eval, wait, php log}. */
const fs = require('fs'), path = require('path'), {execFileSync} = require('child_process');
const {JSDOM, ResourceLoader, VirtualConsole} = require('jsdom');
const BASE = 'http://127.0.0.1:8877/';

function mime(p){ return p.endsWith('.json') ? 'application/json' : p.endsWith('.js') ? 'text/javascript'
  : p.endsWith('.gz') ? 'application/gzip' : p.endsWith('.css') ? 'text/css' : 'text/plain'; }

function makeFetch(root, log){
  return async function fetch(input, init){
    init = init || {};
    const u = new URL(String(input), BASE);
    if (u.origin !== new URL(BASE).origin) return new Response('blocked in harness', {status: 599});
    const rel = decodeURIComponent(u.pathname.replace(/^\//, ''));
    if (rel.endsWith('.php')){
      const get = Object.fromEntries(u.searchParams.entries());
      let post = {};
      if ((init.method || 'GET').toUpperCase() === 'POST' && init.body != null)
        post = Object.fromEntries(new URLSearchParams(String(init.body)).entries());
      const t0 = Date.now();
      const out = execFileSync(process.env.PHP_BIN || 'php', [path.join(__dirname, 'php_shim.php')], {
        env: Object.assign({}, process.env, {SNPT_REQ: JSON.stringify({method: (init.method||'GET').toUpperCase(),
             get, post, root, script: rel})}), maxBuffer: 1 << 28});
      log.push({url: rel, method: init.method || 'GET', ms: Date.now() - t0, bytes: out.length,
                post: Object.keys(post), n_genotypes: post.genotypes ? JSON.parse(post.genotypes).length : undefined,
                reply: out.toString().slice(0, 300)});
      return new Response(out, {status: 200, headers: {'Content-Type': 'application/json'}});
    }
    const f = path.join(root, rel);
    if (!fs.existsSync(f) || fs.statSync(f).isDirectory()) { log.push({url: rel, status: 404}); return new Response('not found', {status: 404}); }
    log.push({url: rel, status: 200});
    return new Response(fs.readFileSync(f), {status: 200, headers: {'Content-Type': mime(f)}});
  };
}

class Loader extends ResourceLoader {
  constructor(root, log){ super(); this.root = root; this.log = log; }
  fetch(url, opts){
    if (/d3js\.org\/d3\.v5/.test(url)) return Promise.resolve(fs.readFileSync(require.resolve('d3/dist/d3.min.js')));
    if (!url.startsWith(BASE)) return Promise.resolve(Buffer.from(''));
    const f = path.join(this.root, decodeURIComponent(new URL(url).pathname.replace(/^\//, '')));
    if (!fs.existsSync(f)){ this.log.push({url, status: 404, kind: 'script'}); return Promise.resolve(Buffer.from('')); }
    return Promise.resolve(fs.readFileSync(f));
  }
}

async function openSite(root){
  const log = [], errors = [];
  const vc = new VirtualConsole();
  vc.on('jsdomError', e => errors.push('jsdomError: ' + (e.message || e)));
  vc.on('error', (...a) => errors.push('console.error: ' + a.map(String).join(' ')));
  vc.on('warn', (...a) => errors.push('console.warn: ' + a.map(String).join(' ')));
  const html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
  const dom = new JSDOM(html, {url: BASE + 'index.html', runScripts: 'dangerously', resources: new Loader(root, log),
    pretendToBeVisual: true, virtualConsole: vc,
    beforeParse(w){
      w.fetch = makeFetch(root, log);
      w.Response = Response; w.Blob = Blob; w.DecompressionStream = DecompressionStream;
      w.TextDecoder = TextDecoder; w.TextEncoder = TextEncoder;
      w.URL.createObjectURL = () => 'blob:harness'; w.URL.revokeObjectURL = () => {};
      w.scrollTo = () => {}; w.HTMLElement.prototype.scrollIntoView = function(){};
      // jsdom has no canvas backend: a no-op 2D context lets the GWAS Explorer's
      // Manhattan plot code run (nothing is drawn; the SNP data and region logic are real).
      const noop = () => {};
      w.HTMLCanvasElement.prototype.getContext = function(){
        const base = {canvas: this, measureText: () => ({width: 0}), getImageData: () => ({data: []}),
                      createLinearGradient: () => ({addColorStop: noop}), getLineDash: () => []};
        return new Proxy(base, {get: (t, k) => (k in t ? t[k] : noop), set: (t, k, v) => { t[k] = v; return true; }});
      };
    }});
  await new Promise(r => dom.window.addEventListener('load', r));
  await new Promise(r => setTimeout(r, 50));
  const $eval = s => dom.window.eval(s);
  const wait = ms => new Promise(r => setTimeout(r, ms));
  return {dom, window: dom.window, $eval, wait, log, errors};
}
module.exports = {openSite};
