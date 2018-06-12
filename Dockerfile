# Use an official g++ runtime as base
FROM gcc:8.1

# Set the working directory to /ReactiveEulerSolver
WORKDIR /ReactiveEulerSolver

# Download CMake and install it
RUN wget https://cmake.org/files/v3.11/cmake-3.11.1-Linux-x86_64.tar.gz \
  && tar xzf cmake-3.11.1-Linux-x86_64.tar.gz \
  && mv cmake-3.11.1-Linux-x86_64/ local/ \
  && rm cmake-3.11.1-Linux-x86_64.tar.gz

# Add the local prefix to the PATH environment variable.
ENV PATH="/ReactiveEulerSolver/local/bin:${PATH}"
ENV LINKER_PATH="/ReactiveEulerSolver/local/lib:${LINKER_PATH}"
ENV LD_LINKER_PATH="/ReactiveEulerSolver/local/lib:${LD_LINKER_PATH}"
ENV CPLUS_INCLUDE_PATH="/ReactiveEulerSolver/local/include:${CPLUS_INCLUDE_PATH}"

# Download and Install Boost, jemalloc and hwloc
RUN apt-get update && apt-get install -y \
  libboost-all-dev \
  libjemalloc-dev \
  libhwloc-dev

# Download and Install HPX
RUN git clone --depth=1 https://github.com/STEllAR-GROUP/hpx.git \
  && mkdir hpx-build \
  && cd hpx-build \
  && cmake ../hpx -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/ReactiveEulerSolver/local -DHPX_WITH_MALLOC=jemalloc -DHPX_WITH_TESTS:BOOL=Off -DHPX_WITH_EXAMPLES:BOOL=Off -DHPX_WITH_DOCUMENTATION:BOOL=Off\
  && make -j2 install \
  && cd ../ \
  && rm -rf hpx*

# Download and Install CGNS
RUN git clone https://github.com/CGNS/CGNS.git \
  && mkdir CGNS-build \
  && cd CGNS-build \
  && cmake ../CGNS -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/ReactiveEulerSolver/local \
  && make install \
  && cd ../ \
  && rm -rf CGNS*

# Download and Configure the ReactiveEulerSolver C++ Code
RUN git clone https://github.com/maikel/ReactiveEulerSolver.git code \
  && mkdir build-release \
  && cd build-release \
  && cmake ../code -DCMAKE_BUILD_TYPE=Release \
  && make
