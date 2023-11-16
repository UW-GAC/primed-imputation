FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
    curl \
    default-jre \
    nodejs

RUN cd /usr/local/bin && \
    curl -sL imputationbot.now.sh | bash

RUN cd /usr/local/bin && \
    curl -L https://github.com/trentm/json/raw/master/lib/json.js > json && \
    chmod 755 json
