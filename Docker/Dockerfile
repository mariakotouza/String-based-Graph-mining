# Install R version 3.5
FROM r-base:3.5.0

# Install Ubuntu packages
RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev/unstable \
    libxt-dev \
    libssl-dev

# Download and install ShinyServer (latest version)
RUN wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb

# Install R packages that are required
# TODO: add further package if you need!
RUN apt-get install r-cran-xml -y 
#RUN R -e "install.packages('XML', repos = 'http://www.omegahat.org/R')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('graph'))"
RUN R -e "install.packages(c('shiny','shinyFiles', 'shinyBS', 'shinyalert', 'shinybusy', 'shinythemes', 'dplyr','stringdist','igraph','visNetwork','CINNA','DT', 'gplots', 'ggraph','graphlayouts','cluster','optrees','shallot','aricode','rfUtilities', 'NMF', 'NMI', 'NetPreProc', 'BiocManager'), repos='http://cran.rstudio.com/')"

# copy the app to the image
RUN mkdir /root/app
COPY app /root/app

COPY Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/app')"]
