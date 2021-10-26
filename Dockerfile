FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
    curl \
    default-jre

RUN cd /usr/local/bin && \
    curl -sL imputationbot.now.sh | bash
