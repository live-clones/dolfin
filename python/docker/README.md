# Dockerfile for running pybind11 tests

To avoid rebulding the DOLFIN C++ interface for every change in the
pybind11 interface, tests are run against a Docker image that is
rebuilt as required.

## Building the image

1. `docker build .`
2. Tags image
  ```
  docker tag my_image fenicsproject/pybind11-testing
  ```
## Uploading to DockerHub

1. `docker login`
2. `docker push fenicsproject/pybind11-testing`
