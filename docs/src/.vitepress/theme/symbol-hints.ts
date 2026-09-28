import nomenclature from '../nomenclature.json'

// ── Symbol hints ────────────────────────────────────────────────────────────
// Hovering an equation lists the symbols of the nomenclature it holds, with
// their meaning on this page: their ids are on the rendered equation
// (`data-nomen`, set when the formula is typeset), and the entries come from the
// same nomenclature.json the typesetting read, with the value of a constant
// read from the library when the documentation was built.

const nomenIndex = new Map((nomenclature as any[]).map((e) => [e.id, e]))

function escapeText(s: string): string {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
}

export function installSymbolHints(): void {
  const tip = document.createElement('div')
  tip.className = 'cl-nomen-hint'
  tip.setAttribute('role', 'tooltip')
  document.body.appendChild(tip)
  let current: Element | null = null

  const hide = () => { tip.classList.remove('is-visible'); current = null }

  document.addEventListener('mouseover', (e) => {
    const eq = (e.target as HTMLElement | null)?.closest?.('mjx-container[data-nomen]')
    if (!eq || eq === current) return
    current = eq
    const rows = (eq.getAttribute('data-nomen') || '').split(' ')
      .map((id) => nomenIndex.get(id)).filter((x) => x)
    if (rows.length === 0) return
    // The label is HTML written in the nomenclature (<sub>, <sup>); the rest is text.
    tip.innerHTML = rows.map((r: any) => {
      const unit = r.value ? `${escapeText(r.value)} ${escapeText(r.unit)}` : (r.unit && r.unit !== '1' ? escapeText(r.unit) : '')
      return `<div class="row"><span class="sym">${r.label}</span><span class="what">${r.name}${unit ? `<span class="unit"> — ${unit}</span>` : ''}</span></div>`
    }).join('')
    tip.classList.add('is-visible')
    const r = eq.getBoundingClientRect()
    const w = Math.min(420, window.innerWidth - 24)
    tip.style.width = `${w}px`
    tip.style.left = `${Math.max(12, Math.min(r.left, window.innerWidth - w - 12))}px`
    tip.style.top = `${r.bottom + window.scrollY + 6}px`
  })

  document.addEventListener('mouseout', (e) => {
    const eq = (e.target as HTMLElement | null)?.closest?.('mjx-container[data-nomen]')
    const to = (e as MouseEvent).relatedTarget as Node | null
    if (eq && eq === current && !(to && eq.contains(to))) hide()
  })
  window.addEventListener('scroll', hide, { passive: true })
}

