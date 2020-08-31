FROM alpine:3.12

RUN apk add --no-cache git gcc libc-dev py3-pip python3-dev

RUN /usr/bin/python3 -m pip install Cython
RUN /usr//bin/python3 -m pip install numpy
RUN /usr/bin/python3 -m pip install git+https://github.com/edwardslab-wustl/dxm.git

