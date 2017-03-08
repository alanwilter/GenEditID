# How to...

## Shiny App

### Installing R 3.3.2
[https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/)

R 3.3.2 binary for Mac OS X 10.9 (Mavericks) and higher, signed package. Contains R 3.3.2 framework, R.app GUI 1.68 in 64-bit for Intel Macs, Tcl/Tk 8.6.0 X11 libraries and Texinfo 5.2. The latter two components are optional and can be ommitted when choosing "custom install", it is only needed if you want to use the tcltk R package or build package documentation from sources.

### Installing RStudio 1.0.136
[https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/)

### How to build a shiny app
- [http://shiny.rstudio.com/tutorial/](http://shiny.rstudio.com/tutorial/)
- [http://shiny.rstudio.com/tutorial/lesson1/](http://shiny.rstudio.com/tutorial/lesson1/)

### How to run the shiny app
- Installation instruction
  - in RStudio, File > New Project and select directory of the git repo
  - install these packages in RStudio
  ```R
  install.packages(c('shiny', 'reshape2', 'ggplot2', 'grofit', 'plotly'))
  ```
- Run the app
```R
source('shinyapp/global.R')
runApp('shinyapp/.')
```
- Choose CSV File, click on 'Browse...' and select all 13 data files in `shinyapp` directory.
