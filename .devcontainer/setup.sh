#!/bin/bash
set -e

# Update package list
apt-get update

# Configure locales
apt-get install -y locales
sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen
locale-gen
echo "export LC_ALL=en_US.UTF-8" >> /etc/bash.bashrc
echo "export LANG=en_US.UTF-8" >> /etc/bash.bashrc
echo "export LANGUAGE=en_US.UTF-8" >> /etc/bash.bashrc

# Install basic tools
apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    wget \
    git \
    build-essential \
    apt-utils \
    software-properties-common \
    lsb-release

# Install Python and common packages
apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-dev \
    python3-venv

# Install pip packages
pip3 install --no-cache-dir \
    numpy \
    pandas \
    matplotlib \
    scipy \
    ipykernel \
    jupyterlab

# Add R repository and key
apt-get install -y --no-install-recommends dirmngr
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

# Install R and common packages
apt-get update
apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev

# Install JAGS
apt-get install -y --no-install-recommends jags

# Set locale for R
echo "LC_ALL=en_US.UTF-8" >> /etc/R/Renviron
echo "LANG=en_US.UTF-8" >> /etc/R/Renviron

# Install R packages
Rscript -e "install.packages(c('rjags', 'coda', 'ggplot2', 'dplyr', 'tidyr', 'IRkernel', 'remotes', 'R2jags'), repos='https://cloud.r-project.org/')"
Rscript -e "IRkernel::installspec(user = FALSE)"

# Install VSCode R debugging tools
Rscript -e "remotes::install_github('ManuelHentschel/vscDebugger', dependencies = TRUE)"

# Create directory for VS Code settings
mkdir -p /home/vscode/.vscode-server/data/Machine/
cat > /home/vscode/.vscode-server/data/Machine/settings.json << EOF
{
    "r.debugger.timeouts.startup": 8000,
    "r.lsp.debug": true,
    "r.rterm.linux": "/usr/bin/R",
    "r.bracketedPaste": true,
    "r.sessionWatcher": true
}
EOF

# Set proper ownership for VS Code settings
chown -R vscode:vscode /home/vscode/.vscode-server

# Clean up
apt-get clean
rm -rf /var/lib/apt/lists/*

echo "Setup complete with R2jags and vscDebugger installation!" 