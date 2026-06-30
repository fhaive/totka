FROM totka_base
COPY renv.lock.prod renv.lock
RUN R -e 'options(renv.config.pak.enabled = FALSE);renv::restore()'
COPY TOTKA_*.tar.gz /app.tar.gz
RUN R -e 'remotes::install_local("/app.tar.gz",upgrade="never")'
RUN tar -xzf /app.tar.gz
RUN rm /app.tar.gz
EXPOSE 80
USER root
RUN sudo chmod -R 777 TOTKA/
USER rstudio
#CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');library(TOTKA);setwd('/TOTKA/');print(getwd());TOTKA::run_app()"
