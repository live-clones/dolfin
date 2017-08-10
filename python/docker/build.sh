#!/bin/bash
docker build --no-cache --tag fenicsproject/pybind11-testing:latest .
docker push fenicsproject/pybind11-testing:latest
