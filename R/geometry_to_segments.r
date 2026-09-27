# Convert Line Geometries to Individual Segments
#
# Extracts linear components from an sf geometry and splits each
# LINESTRING into segments connecting consecutive vertices.
#
# Supports LINESTRING, MULTILINESTRING and mixed GEOMETRYCOLLECTION
# inputs. Empty geometries and zero-length segments are discarded.
# Duplicate and overlapping segments are preserved.
#
# @param geometry An sfc object containing two-dimensional
#   geometries in a planar coordinate system.
#
# @return A numeric matrix with four columns: x0, y0, x1 and y1.
#   Each row represents one segment. If no valid segments exist,
#   an empty matrix with the same columns is returned.
#
# @keywords internal
geometry_to_segments <- function(geometry) {

  # Initialize a consistently structured empty result.
  empty_result <- matrix(
    data = numeric(0),
    ncol = 4L,
    dimnames = list(NULL, c("x0", "y0", "x1", "y1"))
  )

  # Validate the input type.
  if (!inherits(geometry, "sfc")) {
    stop("geometry must be an sfc object.")
  }

  # Return early for an empty input.
  if (length(geometry) == 0L ||
      all(sf::st_is_empty(geometry))) {
    return(empty_result)
  }

  # Require two-dimensional XY geometries.
  dimensions <- vapply(
    geometry,
    function(g) {
      if (sf::st_is_empty(g)) {
        return(TRUE)
      }
      identical(class(g)[1L], "XY")
    },
    logical(1)
  )

  if (!all(dimensions)) {
    stop("Only two-dimensional XY geometries are supported.")
  }

  # Extract linear components only if non-linear geometries are present.
  geometry_types <- sf::st_geometry_type(
    geometry,
    by_geometry = TRUE
  )

  if (!all(geometry_types %in% c("LINESTRING", "MULTILINESTRING"))) {
    geometry <- sf::st_collection_extract(
      geometry, "LINESTRING"
    )
  }

  if (length(geometry) == 0L ||
      all(sf::st_is_empty(geometry))) {
    return(empty_result)
  }

  # Split multipart geometries into individual lines.
  lines <- sf::st_cast(
    geometry, "LINESTRING"
  )

  # Convert consecutive vertices into individual segments.
  segment_list <- lapply(seq_along(lines), function(i) {

    if (sf::st_is_empty(lines[i])) {
      return(NULL)
    }

    xy <- sf::st_coordinates(lines[i])

    if (nrow(xy) < 2L) {
      return(NULL)
    }

    line_segments <- cbind(
      x0 = xy[-nrow(xy), 1L],
      y0 = xy[-nrow(xy), 2L],
      x1 = xy[-1L, 1L],
      y1 = xy[-1L, 2L]
    )

    # Reject invalid coordinates and zero-length segments.
    finite <- rowSums(is.finite(line_segments)) == 4L

    nonzero <- finite &
      (line_segments[, "x0"] != line_segments[, "x1"] |
       line_segments[, "y0"] != line_segments[, "y1"])

    line_segments[nonzero, , drop = FALSE]
  })

  # Combine all segments into a single matrix.
  segment_matrix <- do.call(rbind, segment_list)

  if (is.null(segment_matrix) || nrow(segment_matrix) == 0L) {
    return(empty_result)
  }

  rownames(segment_matrix) <- NULL

  segment_matrix
}
