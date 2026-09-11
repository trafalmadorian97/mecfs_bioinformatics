#import "@preview/cetz:0.3.4"

#set page(width: auto, height: auto, margin: 10pt)
#set text(size: 10pt)

#let teal = rgb("#1b9e77")
#let orange = rgb("#d95f02")
#let ink = rgb("#333333")
#let muted = rgb("#8a8a8a")

#cetz.canvas(length: 1cm, {
  import cetz.draw: *

  // ---- Genealogy (top) ----
  // Present-day haplotypes (tips) and the internal nodes of the tree.
  let H1 = (0, 0); let H2 = (1, 0); let H3 = (2, 0); let H4 = (3, 0); let H5 = (4, 0)
  let L = (0.5, 1.2) // ancestor of H1, H2
  let R2 = (3.5, 0.9) // ancestor of H4, H5
  let R = (2.75, 1.8) // ancestor of H3, {H4,H5}
  let RT = (1.625, 2.8) // root

  let branch(a, b) = line(a, b, stroke: 1pt + ink)
  branch(L, H1); branch(L, H2)
  branch(R2, H4); branch(R2, H5)
  branch(R, H3); branch(R, R2)
  branch(RT, L); branch(RT, R)

  for (p, name) in ((H1, "H1"), (H2, "H2"), (H3, "H3"), (H4, "H4"), (H5, "H5")) {
    content((p.at(0), p.at(1) - 0.3), text(size: 8pt, fill: muted)[#name])
  }

  // A mutation is a coloured dot on the branch where it arose, plus a label.
  let mut(p, col, lbl) = {
    circle(p, radius: 0.12, fill: col, stroke: 0.6pt + white)
    content((p.at(0) + 0.3, p.at(1) + 0.08), text(fill: col, weight: "bold")[#lbl])
  }
  // A and B arise on the SAME branch (root -> ancestor of H1, H2).
  mut((1.175, 2.16), teal, "A")
  mut((0.9, 1.75), teal, "B")
  // C arises on a DIFFERENT branch (ancestor of H3,H4,H5 -> ancestor of H4,H5).
  mut((3.125, 1.35), orange, "C")

  // ---- Sequence: beads on a string (bottom) ----
  // Each bead is a variant at its genomic position; colour marks the branch it
  // arose on. A and B share a colour but sit far apart; C differs from B but is
  // adjacent to it.
  let sy = -1.5
  line((0.2, sy), (4.2, sy), stroke: 1pt + muted)
  let bead(x, col, lbl) = {
    circle((x, sy), radius: 0.18, fill: col, stroke: 0.6pt + white)
    content((x, sy + 0.45), text(fill: col, weight: "bold")[#lbl])
  }
  bead(0.6, teal, "A")
  bead(3.2, orange, "C")
  bead(3.7, teal, "B")

  line((0.2, sy - 0.5), (4.2, sy - 0.5), mark: (end: ">"), stroke: 0.7pt + muted)
  content((2.2, sy - 0.8), text(size: 8pt, style: "italic", fill: muted)[genomic position])

  content((-0.5, 1.4), text(size: 8pt, style: "italic", fill: muted)[genealogy], anchor: "east")
  content((-0.5, sy), text(size: 8pt, style: "italic", fill: muted)[sequence], anchor: "east")
})

#v(6pt)
#block(width: 10.5cm, {
  set text(size: 8.5pt, fill: rgb("#333333"))
  [Variants *A* and *B* arose on the same branch of the genealogy, so the same
    haplotypes carry both: they are perfectly correlated even though they lie far
    apart on the sequence. Variant *C* arose on a different branch and is
    uncorrelated with *A* and *B*, despite sitting right next to *B*. Without
    recombination, correlation reflects shared ancestry, not physical distance.]
})
