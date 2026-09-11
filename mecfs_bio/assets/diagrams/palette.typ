// Shared colour palette for the project's Typst diagrams.
//
// Import the colours you need from a sibling diagram source, e.g.
//   #import "palette.typ": primary, secondary, tertiary, ink, muted
// so that every diagram draws from one consistent palette.
//
// The categorical hues are the Dark2 set, chosen to stay distinguishable for
// the common forms of colour blindness.

#let primary = rgb("#1b9e77") // teal
#let secondary = rgb("#d95f02") // orange
#let tertiary = rgb("#7570b3") // purple
#let quaternary = rgb("#e7298a") // magenta

#let ink = rgb("#333333") // primary line / text colour
#let muted = rgb("#8a8a8a") // secondary lines, axis labels, de-emphasis
