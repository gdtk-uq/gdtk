
#let body-font = "Gelasio"
#let body-size = 9pt
#let heading-font = "Lato"
#let code-font = "Roboto Mono"
#let code-size = 9pt

#let cram-snap(
  title: none,
  icon: none,
  column-number: 2,
  subtitle: none,
  fill-color: "F2F2F2",
  stroke-color: "21222C",
  doc,
) = {

  set text(
    font: body-font,
    size: body-size,
    fill: luma(30),
  )

  show raw: set text(font: code-font, size: code-size)
  show heading: set text(font: heading-font)
  
  let table_stroke(color) = (
    (x, y) => (
      left: none,
      right: none,
      top: none,
      bottom: if y == 0 {
        color
      } else {
        0pt
      },
    )
  )

  let table_fill(color) = (
    (x, y) => {
      if calc.odd(y) {
        rgb(color)
      } else {
        none
      }
    }
  )

  set table(
    align: left + horizon,
    columns: (2.5fr, 3fr),
    fill: table_fill(rgb(fill-color)),
    //stroke: table_stroke(rgb(stroke-color)),
    stroke: none,
  )

  set table.header(repeat: false)

  show table.cell.where(y: 0): set text(weight: "bold", size: 1.2em, font: heading-font)

  columns(column-number)[
    #align(center)[
      #box(height: 1.8em)[
        #if icon != none {
          set image(height: 100%)
          box(icon, baseline: 20%)
          h(0.4cm)
        }
        #text(1.6em, font: heading-font, weight: "bold", tracking: 1.5pt, title)
      ]

      #text(1.6em, font: heading-font, weight: "bold", tracking: 1pt, style: "italic", subtitle)
    ]

    #doc
  ]
}

#let theader(..cells, colspan: 2) = table.header(
  ..cells.pos().map(x => if type(x) == content and x.func() == table.cell {
    x
  } else {
    table.cell(colspan: colspan, x)
  }),
  ..cells.named(),
)
