# Helper function to determine the base directory for a given chapter and section.
# Arguments:
#   cur.dir: The current working directory.
#   sec.dir: The relative path to the chapter and section directory.
#
# Returns:
#   The normalized path to the base directory for the specified chapter and section.
# Example:
#   get.base.dir(getwd(), "Chapter 5/Sec 5.2.1/")

get.base.dir <- function(cur.dir, sec.dir) {
  if (!is.character(cur.dir) || length(cur.dir) != 1L || is.na(cur.dir)) {
    stop("cur.dir must be one non-missing string.")
  }

  if (!is.character(sec.dir) || length(sec.dir) != 1L || is.na(sec.dir)) {
    stop("sec.dir must be one non-missing string.")
  }

  cur.dir <- normalizePath(cur.dir, winslash = "/", mustWork = TRUE)
  sec.dir <- gsub("\\\\", "/", sec.dir)
  sec.parts <- strsplit(sec.dir, "/", fixed = TRUE)[[1]]
  sec.parts <- sec.parts[nzchar(sec.parts)]

  if (length(sec.parts) == 0L) {
    stop("sec.dir must identify a chapter and section directory.")
  }

  cur.parts <- strsplit(cur.dir, "/", fixed = TRUE)[[1]]
  already.in.section <- length(cur.parts) >= length(sec.parts) &&
    identical(
      tolower(tail(cur.parts, length(sec.parts))),
      tolower(sec.parts)
    )

  if (already.in.section) {
    base.dir <- cur.dir
  } else {
    base.dir <- file.path(cur.dir, sec.dir)
  }

  normalizePath(base.dir, winslash = "/", mustWork = TRUE)
}
