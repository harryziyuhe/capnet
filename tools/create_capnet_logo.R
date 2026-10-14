# Reproducible capnet R package hexagon logo (concept 3.1)
# Uses grid plus ragg for reliable PNG output. svglite is optional for SVG.
# Run from the package root: Rscript tools/create_capnet_logo.R

capnet_logo <- function(
  png_file = "man/figures/logo.png",
  svg_file = "man/figures/logo.svg",
  size_px = 1200L
) {
  stopifnot(is.character(png_file), length(png_file) == 1L,
            is.character(svg_file), length(svg_file) == 1L,
            is.numeric(size_px), length(size_px) == 1L,
            is.finite(size_px), size_px >= 300)

  # Fixed, platform-independent geometry in normalized [0, 1] coordinates.
  # No random sampling is used; rerunning produces the same design.
  navy      <- "#102944"
  blue      <- "#4888C0"
  red       <- "#E94D46"
  pale_blue <- "#A9C0D0"
  gray_blue <- "#C8D4DF"
  white     <- "#FFFFFF"

  draw <- function() {
    grid::grid.newpage()

    # Point-up, classic R-package hexagon.
    hx <- c(0.50, 0.90, 0.90, 0.50, 0.10, 0.10)
    hy <- c(0.97, 0.76, 0.24, 0.03, 0.24, 0.76)
    grid::grid.polygon(
      x = hx, y = hy, default.units = "npc",
      gp = grid::gpar(fill = white, col = navy, lwd = 6.5,
                      linejoin = "round")
    )

    # Pale trajectories stay near the fitted path but exceed both limits.
    t <- seq(0, 1, length.out = 300)
    x <- 0.155 + 0.690 * t
    envelope <- sin(pi * t)
    y_controlled <- 0.49 + 0.13 * t + 0.060 * sin(3 * pi * t)
    gray_curves <- list(
      y_controlled + 0.24 * envelope * sin(2 * pi * t + 0.1 * pi),
      y_controlled + 0.22 * envelope * sin(2 * pi * t + 0.9 * pi),
      y_controlled + 0.20 * envelope * sin(2 * pi * t + pi)
    )

    for (i in seq_along(gray_curves)) {
      grid::grid.lines(
        x = x, y = gray_curves[[i]], default.units = "npc",
        gp = grid::gpar(
          col = grDevices::adjustcolor(
            if (i == 2L) pale_blue else gray_blue,
            alpha.f = 0.72
          ),
          lwd = 4.4, lineend = "round", linejoin = "round"
        )
      )
    }

    # Stylized lower and upper contribution thresholds (soft, not hard caps).
    for (y_cap in c(0.405, 0.765)) {
      grid::grid.lines(
        x = c(0.15, 0.85), y = c(y_cap, y_cap),
        default.units = "npc",
        gp = grid::gpar(col = red, lwd = 3.2, lty = "dashed",
                        lineend = "round")
      )
    }

    # Emphasized, stabilized contribution trajectory.
    grid::grid.lines(
      x = x, y = y_controlled, default.units = "npc",
      gp = grid::gpar(col = navy, lwd = 5.5,
                      lineend = "round", linejoin = "round")
    )

    # Two-color logotype, raised to sit closer to the bounded signal. Use a consistent installed font if exact
    # cross-machine text layout is essential; 'sans' is portable.
    grid::grid.text(
      "cap", x = 0.50, y = 0.275, default.units = "npc",
      just = c("right", "centre"),
      gp = grid::gpar(col = navy, fontfamily = "sans",
                      fontface = "bold", fontsize = 65)
    )
    grid::grid.text(
      "net", x = 0.50, y = 0.275, default.units = "npc",
      just = c("left", "centre"),
      gp = grid::gpar(col = blue, fontfamily = "sans",
                      fontface = "bold", fontsize = 65)
    )
  }

  dir.create(dirname(svg_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(png_file), recursive = TRUE, showWarnings = FALSE)
  outputs <- character()

  # SVG: scalable source asset when svglite is installed.
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = svg_file, width = 6, height = 6,
                     bg = "transparent")
    tryCatch(draw(), finally = grDevices::dev.off())
    outputs <- c(outputs, svg = svg_file)
  } else {
    message("Skipping SVG output: install the optional 'svglite' package.")
  }

  # PNG: use ragg when available, avoiding a system X11/Cairo dependency.
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = png_file, width = size_px, height = size_px,
                  units = "px", background = "transparent",
                  res = size_px / 6)
  } else {
    grDevices::png(filename = png_file, width = size_px, height = size_px,
                   bg = "transparent", res = size_px / 6)
  }
  tryCatch(draw(), finally = grDevices::dev.off())
  outputs <- c(outputs, png = png_file)

  invisible(outputs)
}

# Run this script once (or whenever you need to regenerate the assets).
capnet_logo()
