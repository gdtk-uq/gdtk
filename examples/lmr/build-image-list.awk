# RJG, 2025-01-02
#
# NOTE: SVG_ONLY variable can be set to true
#       on command line using: -v SVG_ONLY=1
#
BEGIN {
  imagesDir = ""
}
$1 ~ /:imagesdir:/ {
  imagesDir = $2
}
$1 ~ /^image::/ {
  split($1, a, /::|\[/)
  file = a[2]
  split(file, a, /\./)
  ext = a[2]
  if (SVG_ONLY) {
    if (ext == "svg") 
      print imagesDir"/"file
  }
  else {
    if (ext != "svg")
       print imagesDir"/"file
  }
}
