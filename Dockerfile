FROM ubuntu:22.04

RUN apt-get -qq update \
  && DEBIAN_FRONTEND=noninteractive apt-get install -yq \
  python3-dev \
  python3-pip \
  && \
  rm -rf /var/lib/apt/lists/*

ADD . /opt/laytr-source
WORKDIR /opt/laytr-source

RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install setproctitle pylint anybadge coverage && \
    python3 -m pip install --upgrade setuptools && \
    python3 -m pip install ./

WORKDIR /data

ENTRYPOINT ["laytr"]
