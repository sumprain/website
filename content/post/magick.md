+++
title = "Batch modification of PDF pages with magick package"
date = 2018-01-06T14:46:52+05:30
draft = false

# Tags and categories
# For example, use `tags = []` for no tags, or the form `tags = ["A Tag", "Another Tag"]` for one or more tags.
tags = ["R", "PDF", "magick", "imagemagick"]
categories = ["R"]

# Featured image
# Place your image in the `static/img/` folder and reference its filename below, e.g. `image = "example.jpg"`.
# Use `caption` to display an image caption.
#   Markdown linking is allowed, e.g. `caption = "[Image credit](http://example.org)"`.
# Set `preview` to `false` to disable the thumbnail in listings.
[header]
image = ""
caption = ""
preview = true

+++

## Introduction

I was recently working on a PDF file with 60 odd pages and 
I had to crop the lower border of the pages.

## ImageMagick

[ImageMagick](https://www.imagemagick.org/) is GNU image 
editing program in Linux. Recently a package named 
[magick](https://github.com/ropensci/magick) has been created 
by [rOpenSci](https://github.com/ropensci).

## Steps of batch manipulation of PDF pages

### Loading PDF pages

The PDF pages are loaded as vector of image object

```r
pdf <- image_read("path_file.pdf")
```

We can access page `i` with

```r
pdf[i]
```

### Clipping each page as batch

We run the following command to crop each page from below

```r
pdf <- image_chop(pdf, geometry = "See documentation for the notation")
```

### Save as the new file

```r
image_write(pdf, "path")
```

## Details of the magick package and resources

The documentation for magick is available [here](https://cran.r-project.org/web/packages/magick/vignettes/intro.html), 
[here](https://ropensci.org/blog/2016/08/23/z-magick-release/)

The documentation for ImageMagick is available [here](http://www.imagemagick.org/Usage/).

Comments are welcome.


