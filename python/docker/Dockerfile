# Builds a Docker image with FEniCS stable version built from
# git sources. The image is at:
#
#    https://quay.io/repository/fenicsproject/stable
#
# Authors:
# Chris Richardson

FROM quay.io/fenicsproject/dev-env:latest
MAINTAINER fenics-project <fenics-support@googlegroups.com>

USER root

# Python2 build.
ENV FENICS_BUILD_TYPE=Developer \
    FENICS_PREFIX=/usr/local \
    FENICS_VERSION=master \
    CMAKE_EXTRA_ARGS="-DDOLFIN_ENABLE_PYTHON=False"

RUN /bin/bash -c "FENICS_SRC_DIR=/tmp/src /home/fenics/bin/fenics-pull && cd /tmp/src/dolfin && git checkout garth/feature-pybind11 && git log  | head"

COPY bin $FENICS_HOME/bin

# Install fenics as root user into /usr/local then remove the fenics-* scripts
# the fenics.env.conf file and the unnecessary /home/fenics/local directory as
# the user does not need them in the stable image!
RUN /bin/bash -c "FENICS_SRC_DIR=/tmp/src /home/fenics/bin/fenics-build && \
               ldconfig && \
               rm -rf /home/fenics/local && \
               rm -rf /tmp/src && \
               rm -rf $FENICS_HOME/bin && \
               echo '' >> $FENICS_HOME/.profile"

RUN rm -rvf /usr/local/lib/python3 && sudo pip3 install pkgconfig pytest six --upgrade

