# syntax=docker/dockerfile:1
# MULTI-STAGE DOCKERFILE for building an R package image
# STAGE 1: base_stage is the base image for CI and production
# STAGE 2: development_only_stage is optimised for VSCode development

########################
# STAGE 1: base_stage  #
########################

# IMPORTANT
# If you change the base image, you will need to update the
# PRE_FETCH_BASE_IMAGE variable in the .gitlab-ci.yml file also.
FROM rocker/r-ver:4.4.0 AS base_stage

USER root
# Arguments for UID and GID matching
ARG USER_ID=1000
ARG GROUP_ID=1000

# Set the top level environment variables
ENV \
    DATA_DIRECTORY="/data" \
    OPT_DIRECTORY="/opt" \
    USER_NAME="admin" \
    USER_DIRECTORY="/home/admin" \
    LC_ALL="en_US.UTF-8" \
    LANG="en_US.UTF-8" \
    PKGTYPE="binary"

# Set next environment variables
ENV \
    USER_BASHRC="${USER_DIRECTORY:?}/.bashrc" \
    USER_BIN_DIRECTORY="${USER_DIRECTORY:?}/.local/bin" \
    PROJECT_DIRECTORY="${OPT_DIRECTORY:?}/repo"

# Create user and directories
RUN \
    locale-gen "${LANG:?}" \
    && update-locale LANG="${LANG:?}" \
    && groupadd -g ${GROUP_ID} "${USER_NAME}" \
    && useradd -u ${USER_ID} -g ${GROUP_ID} "${USER_NAME}" --shell /bin/bash --create-home --home-dir "${USER_DIRECTORY}" \
    && if ! getent group docker > /dev/null; then groupadd docker; fi \
    && usermod -a -G docker,staff admin \
    && mkdir -p "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${OPT_DIRECTORY:?}" \
    && chown -R "${USER_NAME:?}:${USER_NAME:?}" "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}" \
    && chmod -R 755 "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}"

################################################
# Install the system dependencies              #
################################################
RUN \
    apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
    build-essential \
    curl \
    wget \
    vim \
    git \
    qpdf \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    texlive-latex-base \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-latex-extra \
    lmodern \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

####################################################
# Install R packages                               #
####################################################

WORKDIR $PROJECT_DIRECTORY

# Copy package files
COPY --chown="${USER_NAME}:${USER_NAME}" ["DESCRIPTION", "NAMESPACE", "./"]

# Install BiocManager and Bioconductor dependencies
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" \
    && R -e "BiocManager::install(c('DESeq2', 'DEGreport', 'SummarizedExperiment'), ask = FALSE, update = FALSE)"

# Install CRAN dependencies and dev tools
RUN R -e "install.packages(c('devtools', 'pkgdown', 'GGally', 'dplyr', 'ggplot2', 'purrr', 'tibble', 'tidyr', 'rlang'), repos='https://cloud.r-project.org')"

# Copy the rest of the project files
COPY --chown="${USER_NAME}:${USER_NAME}" . .

# Install the package
RUN R -e "devtools::install(dependencies = TRUE)"

# Reapply permissions
RUN \
    chown -R "${USER_NAME:?}:${USER_NAME:?}" "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}" \
    && chmod -R 755 "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}"

USER "${USER_NAME:?}"
WORKDIR ${PROJECT_DIRECTORY}

###################################
# STAGE 2: development_only_stage #
###################################

FROM base_stage AS development_only_stage
USER root

# Optional sudo for development
ARG HAS_SUDO="${HAS_SUDO:-0}"
RUN if [ "${HAS_SUDO}" = "1" ]; then \
    apt-get update -y \
    && apt-get install -y sudo \
    && rm -rf /var/lib/apt/lists/* \
    && echo "${USER_NAME:?} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers; \
    fi

WORKDIR "${PROJECT_DIRECTORY}"
USER "${USER_NAME}"

# Prepare VSCode server directories
RUN mkdir -p "${USER_DIRECTORY}/.vscode-server/extensions" \
    "${USER_DIRECTORY}/.vscode-server-insiders/extensions" \
    && chown -R "${USER_NAME}:${USER_NAME}" \
    "${USER_DIRECTORY}/.vscode-server" \
    "${USER_DIRECTORY}/.vscode-server-insiders"
