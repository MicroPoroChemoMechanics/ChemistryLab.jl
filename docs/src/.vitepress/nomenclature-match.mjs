// Which symbols of the nomenclature a formula holds, read from its TeX source.
//
// The entries come from `nomenclature.json`, which make.jl writes from
// docs/nomenclature.toml. The plugin that typesets the formulas calls
// `symbolsOf(tex, page)` on each one and stores the ids it returns on the
// rendered equation; the theme shows those entries when the equation is hovered.
// Plain JavaScript, so that node can test it without the build.

// Commands whose argument is text or a chemical formula, not symbols.
const DROPPED = new Set([
  'ce', 'text', 'textrm', 'textbf', 'textit', 'mathrm', 'operatorname', 'mbox', 'label', 'tag',
])
// Commands that make one symbol of their argument: \dot{B}, \mathcal{A}.
const ACCENTS = new Set([
  'mathcal', 'mathbf', 'boldsymbol', 'mathit', 'mathsf', 'mathbb', 'dot', 'ddot', 'bar',
  'hat', 'tilde', 'vec', 'overline', 'mathring', 'widehat', 'widetilde',
])

// The next atom of `s` from `i`: a braced group (its content), a command, or a
// character. Returns [text, next index].
function atom(s, i) {
  while (i < s.length && /\s/.test(s[i])) i++
  if (i >= s.length) return ['', i]
  if (s[i] === '{') {
    let depth = 0
    let j = i
    for (; j < s.length; j++) {
      if (s[j] === '{') depth++
      else if (s[j] === '}') { depth--; if (depth === 0) break }
    }
    return [s.slice(i + 1, j), j + 1]
  }
  if (s[i] === '\\') {
    const m = /^\\([A-Za-z]+|.)/.exec(s.slice(i))
    return [m[0], i + m[0].length]
  }
  return [s[i], i + 1]
}

/**
 * The symbols of a TeX string, in order: `{ t, script }`, `t` a command
 * (`\sigma`), an accented symbol (`\dot{B}`) or a character, `script` true
 * inside a subscript or superscript.
 */
export function tokenize(tex, inScript = false, out = []) {
  let i = 0
  while (i < tex.length) {
    const c = tex[i]
    if (/\s/.test(c) || c === '{' || c === '}' || c === '&' || c === "'") { i++; continue }
    if (c === '_' || c === '^') {
      const [a, j] = atom(tex, i + 1)
      tokenize(a, true, out)
      i = j
      continue
    }
    if (c === '\\') {
      const m = /^\\([A-Za-z]+)/.exec(tex.slice(i))
      if (!m) { i += 2; continue } // \, \; \! \{ and the like
      const name = m[1]
      i += m[0].length
      if (DROPPED.has(name)) { const [, j] = atom(tex, i); i = j; continue }
      if (ACCENTS.has(name)) {
        const [a, j] = atom(tex, i)
        out.push({ t: `\\${name}{${a.replace(/\s+/g, '')}}`, script: inScript })
        i = j
        continue
      }
      out.push({ t: `\\${name}`, script: inScript })
      continue
    }
    out.push({ t: c, script: inScript })
    i++
  }
  return out
}

function holdsOn(entry, page) {
  return entry.pages.some((p) => (p.endsWith('/') ? page.startsWith(p) : page === p))
}

/**
 * The entries meant on `page` (a path under docs/src, such as
 * "theory/surface_complexation.md"): for each TeX form, the entries scoped to
 * the page if there are any, the unscoped one otherwise. A subscript and a
 * symbol written alike are two forms: the `r` of Δ_r is not the rate `r`.
 */
export function entriesFor(entries, page) {
  const byTex = new Map()
  for (const e of entries) {
    const key = `${e.script ? 'script:' : ''}${e.tex}`
    if (!byTex.has(key)) byTex.set(key, [])
    byTex.get(key).push(e)
  }
  const chosen = []
  for (const group of byTex.values()) {
    const scoped = group.filter((e) => e.pages.length > 0 && holdsOn(e, page))
    chosen.push(...(scoped.length > 0 ? scoped : group.filter((e) => e.pages.length === 0)))
  }
  return chosen
}

/**
 * The ids of the entries a formula holds, in the order they first appear. A
 * longer symbol is matched first and its characters are then taken: the `A` of
 * `N_A` is not also the conservation matrix. An entry marked `script` matches
 * only in a subscript or superscript, any other one only outside.
 */
export function symbolsOf(tex, page, entries) {
  const toks = tokenize(tex)
  const taken = new Array(toks.length).fill(false)
  const found = []
  const candidates = entriesFor(entries, page)
    .map((e) => ({ e, seq: tokenize(e.tex) }))
    .filter(({ seq }) => seq.length > 0)
    .sort((a, b) => b.seq.length - a.seq.length)
  for (const { e, seq } of candidates) {
    for (let p = 0; p + seq.length <= toks.length; p++) {
      let ok = true
      for (let k = 0; k < seq.length && ok; k++) {
        const tok = toks[p + k]
        const script = k === 0 && e.script ? true : seq[k].script
        ok = !taken[p + k] && tok.t === seq[k].t && tok.script === script
      }
      if (!ok) continue
      for (let k = 0; k < seq.length; k++) taken[p + k] = true
      found.push({ id: e.id, at: p })
    }
  }
  const seen = new Set()
  return found
    .sort((a, b) => a.at - b.at)
    .map((f) => f.id)
    .filter((id) => (seen.has(id) ? false : (seen.add(id), true)))
}
