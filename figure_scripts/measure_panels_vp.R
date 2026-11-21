# 视口法测量 patchwork 面板 bbox（绕开本机 null 单位换算失效问题）
# 供 fig 脚本 source 使用：measure_patchwork_panels(plot, width_in, height_in, manifest_path)
measure_patchwork_panels <- function(plot, width_in, height_in, manifest_path) {
  g <- patchwork::patchworkGrob(plot)
  probe <- tempfile(fileext = ".pdf")
  grDevices::pdf(probe, width = width_in, height = height_in, useDingbats = FALSE)
  grid::grid.newpage(); grid::grid.draw(g); grid::grid.force()
  vps <- grid::grid.ls(viewports = TRUE, grobs = FALSE)$name
  pn <- vps[grepl("^panel-[0-9]+[.]", vps)]
  panels <- lapply(pn, function(nm) {
    grid::downViewport(nm)
    bl <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc"))
    tr <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(1, "npc"))
    grid::upViewport(0)
    list(name = nm, bbox = c(grid::convertX(bl$x, "pt", valueOnly = TRUE),
                             grid::convertY(bl$y, "pt", valueOnly = TRUE),
                             grid::convertX(tr$x, "pt", valueOnly = TRUE),
                             grid::convertY(tr$y, "pt", valueOnly = TRUE)))
  })
  grDevices::dev.off(); unlink(probe)
  ids <- letters[seq_along(panels)]
  manifest <- list(
    schema_version = 1L,
    backend = "r-patchwork",
    figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
    panels = lapply(seq_along(panels), function(i) {
      list(id = ids[i], bbox_pt = as.numeric(panels[[i]]$bbox),
           grid_id = "patchwork-grid-1",
           row_start = 0L, row_stop = 1L,
           col_start = as.integer(i - 1), col_stop = as.integer(i))
    }),
    exemptions = list()
  )
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, digits = 6, pretty = TRUE)
  invisible(manifest)
}
