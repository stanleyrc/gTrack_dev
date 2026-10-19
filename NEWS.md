gTrack 0.1.0.9000 (feature/pub-style)
----------------------------------------------------------

* `pub_style()` / `plot(gt, win, style = pub_style())` / `options(gTrack.style = pub_style())`: publication style
  (font family, absolute pt sizes, thin axes with outward ticks of fixed mm length, no grey background or y grid,
  one clean Mb x axis with the chromosome name, left-aligned track titles with optional `name.col`, edge labels
  coloured like their edge at the arc apex, margins in mm, multi-window plots without the empty space reserved per window).
* `gtrack_grob()`: render a gTrack plot as a vector grid grob (gridGraphics) for grid / patchwork layouts.
* `font.label` formatting field (eg 3 = italic gene names).
* Without a style gTrack draws exactly as before.

gTrack 0.1.0 2016-3-18
----------------------------------------------------------

* Initial stable release to Github
