gTrack 0.1.0.9000 (feature/pub-style)
----------------------------------------------------------

* `pub_style()` / `plot(gt, win, style = pub_style())` / `options(gTrack.style = pub_style())`: publication style
  (font family, absolute pt sizes, thin axes with outward ticks of fixed mm length, no grey background or y grid,
  one clean Mb x axis with the chromosome name, left-aligned track titles with optional `name.col`, edge labels
  coloured like their edge at the arc apex, margins in mm, multi-window plots without the empty space reserved per window).
* `gtrack_grob()`: render a gTrack plot as a vector grid grob (gridGraphics) for grid / patchwork layouts.
* `font.label` formatting field (eg 3 = italic gene names).
* `pub_style(win.widths = )`: fixed display share per window (data stay linear inside a window); `xaxis.chr.merge`: one chromosome label per run of windows.
* `pub_style(arc.style = "bell")`: one renderer for every junction in a multi-window plot: ALT edges are bells from both segment ends with height ~ log10(`edges$span`),
  kept inside the track; REF edges straight; edges to segments outside the windows are short arrows labelled `edges$dangle.label` in a row under the track;
  edge labels above the apex with a white halo, moved off other labels and arcs; walk paths on one line joined by bells.
* `pub_style(arc.strand = TRUE, arc.hook.mm = )`: strand-aware bells: each end leaves its breakend sideways towards the joined segment end (dup-like hooks outward, deletion-like arches between inner edges, inversions same side); walk-path bells likewise.
* `name.pos = "right"` formatting: track title right of the windows. `pub_window_layout()`: window positions (fraction of the plot width) for overlays / hit regions.
* Without a style gTrack draws exactly as before.

gTrack 0.1.0 2016-3-18
----------------------------------------------------------

* Initial stable release to Github
