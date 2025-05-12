FROM rocker/r2u:jammy

LABEL org.opencontainers.image.authors="Thomas De Almeida <thomas.de_almeida@ens-lyon.fr>, Samuel Dunesme <samuel.dunesme@ens-lyon.fr>"
LABEL org.opencontainers.image.source="https://github.com/tdealmeida/RiverSeg"
LABEL org.opencontainers.image.documentation="https://tdealmeida.github.io/RiverSeg/"
LABEL org.opencontainers.image.description="Shiny application to visualize segmentation to visualize results both graphically and cartographically according to a set of segmentation methods "

RUN locale-gen fr_FR.UTF-8

RUN Rscript -e 'install.packages("shiny")'
RUN Rscript -e 'install.packages("shinyBS")'
RUN Rscript -e 'install.packages("changepoint")'
RUN Rscript -e 'install.packages("dplyr")'
RUN Rscript -e 'install.packages("RPostgres")'
RUN Rscript -e 'install.packages("ecp")'
RUN Rscript -e 'install.packages("cpm")'
RUN Rscript -e 'install.packages("reticulate")'
RUN Rscript -e 'install.packages("zoo")'
RUN Rscript -e 'install.packages("Rbeast")'
RUN Rscript -e 'install.packages("readr")'
RUN Rscript -e 'install.packages("sf")'
RUN Rscript -e 'install.packages("leaflet")'
RUN Rscript -e 'install.packages("ggplot2")'
RUN Rscript -e 'install.packages("EnvCpt")'
RUN Rscript -e 'install.packages("kcpRS")'
RUN Rscript -e 'install.packages("cumSeg")'
RUN Rscript -e 'install.packages("trend")'
RUN Rscript -e 'install.packages("bfast")'

RUN Rscript -e 'install.packages("remotes")'
RUN R -e 'remotes::install_github("lvaudor/hubr")'

RUN mkdir /app
ADD . /app
WORKDIR /app

EXPOSE 3840

RUN groupadd -g 1010 app && useradd -c 'app' -u 1010 -g 1010 -m -d /home/app -s /sbin/nologin app
USER app

CMD  ["R", "-e", "shiny::runApp('.', port=3840, host='0.0.0.0')"]