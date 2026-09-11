#import "@preview/cetz:0.3.4"
#import "palette.typ": ink, muted, primary, secondary, tertiary

#set page(width: auto, height: auto, margin: 10pt)
#set text(size: 10pt)

#cetz.canvas(length: 1cm, {
  import cetz.draw: *

  // ---- Genealogy (left panel) ----
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
  // Variants sharing a branch share a colour (they are perfectly correlated).
  let mut(p, col, lbl, dx: 0.28, dy: 0.06) = {
    circle(p, radius: 0.11, fill: col, stroke: 0.6pt + white)
    content((p.at(0) + dx, p.at(1) + dy), text(size: 9pt, fill: col, weight: "bold")[#lbl])
  }
  // A, B on the branch root -> ancestor of H1, H2.
  mut((1.175, 2.16), primary, "A")
  mut((0.9275, 1.808), primary, "B")
  // C, D on the branch ancestor(H3,H4,H5) -> ancestor of H4, H5.
  mut((3.05, 1.44), secondary, "C")
  mut((3.2375, 1.215), secondary, "D")
  // E, F, G on the branch leading to H3.
  mut((2.525, 1.26), tertiary, "E", dx: -0.3)
  mut((2.375, 0.9), tertiary, "F", dx: -0.3)
  mut((2.225, 0.54), tertiary, "G", dx: -0.3)

  // ---- Sequence: beads on a string (right panel) ----
  // Each bead is a variant at its genomic position; colour marks correlation
  // group. Same-colour beads sit at varying distances -- distance is irrelevant.
  let sy = 1.4
  let (sx0, sx1) = (5.8, 10.7)
  line((sx0, sy), (sx1, sy), stroke: 1pt + muted)
  let bead(x, col, lbl) = {
    circle((x, sy), radius: 0.15, fill: col, stroke: 0.6pt + white)
    content((x, sy + 0.42), text(size: 9pt, fill: col, weight: "bold")[#lbl])
  }
  bead(5.95, primary, "A")
  bead(6.65, tertiary, "E")
  bead(7.45, secondary, "C")
  bead(8.2, tertiary, "F")
  bead(9.0, primary, "B")
  bead(9.7, secondary, "D")
  bead(10.5, tertiary, "G")

  // Panel titles.
  content((2.0, 3.3), text(size: 9pt, style: "italic", fill: muted)[genealogy])
  content((8.25, 3.3), text(size: 9pt, style: "italic", fill: muted)[sequence])
})

#v(6pt)
#align(center, {
  set text(size: 9pt, fill: ink)
  [Dots indicate variants. Dots with the same colour are perfectly correlated ($r^2 = 1$).]
})
