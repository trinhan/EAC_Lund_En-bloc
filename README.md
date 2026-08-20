# EAC Lund En-bloc

This repository contains public-facing analysis code and final interpreted result summaries for the EAC Lund en-bloc cohort analysis.

The bookdown link can be found in https://trinhan.github.io/EAC_Lund_En-bloc/

Raw sequencing files, raw variant calls, clinical metadata, and private cohort governance materials are not included. The public tables in `public_outputs/` are derived from the private analysis workflow and are limited to summarized genomic results intended for display.

The Bookdown site can be rendered from the sanitized files in `public_outputs/`; rendered GitBook files are written to `book/`.

```sh
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"
```


