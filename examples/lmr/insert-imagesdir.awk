# RJG, 2025-01-01
#
# This script is used as part of a pipeline to construct a top-level listing
# of examples in lmr. It is used to transform the output from asciidoctor-reducer.
# asciidoctor-reducer is first used to collate the "description" portions
# from READMEs on a per example basis. This AWK script can then set the
# :imagesdir: appropriate for each example. In this way, image files are
# located correctly when building the HTML.

BEGIN {
  PAT_TO_MATCH = "`gdtk/examples/lmr/"
  LEN_PAT = length(PAT_TO_MATCH)
}
$0 ~ PAT_TO_MATCH {
  imagesdir = substr($1, LEN_PAT+1)
  imagesdir = substr(imagesdir, 1, length(imagesdir)-1)
  print $0
  print ""
  print ":imagesdir: "imagesdir
  next
}
1 {
  print $0
}
