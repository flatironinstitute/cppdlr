# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang
ARG APPNAME=app4triqs

# Install here missing dependencies, e.g.
RUN apt-get update && apt-get install -y libzstd-dev

COPY --chown=build . $SRC/$APPNAME
RUN mkdir $BUILD/$APPNAME && chown build $BUILD/$APPNAME

RUN if echo "$CXX" | grep -q "clang"; then \
        git clone https://github.com/flatironinstitute/clair --branch unstable $SRC/clair && \
        mkdir $BUILD/clair && \
        CXXFLAGS="" cmake -S $SRC/clair -B $BUILD/clair -DCMAKE_INSTALL_PREFIX=$INSTALL $CMAKE_ARGS && \
        cmake --build $BUILD/clair && cmake --install $BUILD/clair; \
    fi

ARG BUILD_ID
ARG CMAKE_ARGS
USER build
WORKDIR $BUILD/$APPNAME
RUN cmake $SRC/$APPNAME -DTRIQS_ROOT=${INSTALL} $CMAKE_ARGS && make -j4 || make -j1 VERBOSE=1
USER root
RUN make install
