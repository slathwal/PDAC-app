# Introduction
This repository contains the source code for R Shiny application hosted at [slathwal.shinyapps.io/pdac-app](https://slathwal.shinyapps.io/pdac-app/). The app contains analysis of pancreatic cancer data from the Cancer Genome Atlas (TCGA).

## Why was the application created?
The application was created to play with some of the new R packages and to keep my skills up-to-date while on a break. Some interesting packages used in the app are as follows:
- [rtables](https://insightsengineering.github.io/rtables/main/index.html) - A package developed as part of work done by working groups at [R consortium](https://r-consortium.org/). The package allows creation of tables that are ready for submission of clinical trial data to regulatory health authorities.
- [circlize](https://jokergoo.github.io/circlize/) - A package used to create data visualizations in a circular format. The package is very complex and I used only a part of it. ChatGPT was instrumental in helping me understand and get to the parts I needed quickly.
- [shinywidgets](https://dreamrs.github.io/shinyWidgets/index.html) - A versatile package that I used to display informational notes to users.

## Additional resources
The following resources were very helpful while creating the application:
- [Mastering Shiny](https://mastering-shiny.org/index.html)
- [Advanced R](https://adv-r.hadley.nz/)