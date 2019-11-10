# SGM: String-based Graph Mining Framework

SGM is a software framework that provides visual analytics services on string-based data, implemented in [R Shiny](https://shiny.rstudio.com/).
It is an explorative and interactive framework that provides graph construction, graph sparsification, graph management 
and graph mining processes, including central node extraction, node clustering, classification and evaluation.

## Installation

#### R Shiny packages

```
install.packages(c("shiny","shinyFiles", "shinyBS", "shinyalert", "shinybusy", "shinythemes"))
```

#### Other packages for data processing and visualization

```
install.packages(c("dplyr","stringdist","igraph","visNetwork","CINNA","DT", "gplots", "ggraph","graphlayouts","cluster","optrees","shallot","aricode","rfUtilities", "NMF", "NMI", "NetPreProc", "graph"))
```

##  Run the Application

- You can run the R shiny application (**Shiny-App** folder) by pressing the **"Run App"** button, from either **ui.R** or **server.R** script.
- Before running the app, you may change domain variable in **global.R** script, according to the application domain. For gene-sequence data for example, 'AA' is used as domain name.
- In case that your data are not in script-based format, you may use the data transformation module to transform them.
- Sample datasets from two application domains: 1) gene-sequence data and 2) document data , are given in the **Dataset** folder (coming soon)

##  Run SGM as an R script-based tool 
In order to run SGM as a script-based tool, except from the packages that are described above, the user need to install the optparse package.

```
install.packages("optparse") 
```

There are two ways to run the script-based tool:
- Through R Studio: run the **make_options.R** file, after first changing/editing the default values for the parameters that are in the **option_list** in the file. The working directory should be the path where the file **make_options.R** is.
- Through the command line: run the command **Rscript --vanilla make_options.R** followed by a list of the parameters you need to change from the default values. For example, in order to run only the 1st & 2nd pipeline choice, the command should look like this:

```
Rscript --vanilla make_options.R --pipeline 1,2 
```

All the available and the deafault parameters are available in the **make_options.R** files. The user can ask for help using the following command:
```
Rscript --vanilla make_options.R --help
```

##  Run SGM as a Docker container
The Dockerfile and the configurations needed for building the SGM docker image are available in the *docker* folder.

The docker image of SGM is available on DockerHub throught the following link:
https://cloud.docker.com/u/mariakotouza/repository/docker/mariakotouza/sgm.

##  License
This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International  [LICENSE](License.md). For more details visit https://creativecommons.org/licenses/by-nc-sa/4.0/.