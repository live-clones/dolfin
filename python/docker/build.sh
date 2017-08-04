#!/bin/bash
docker build --tag fenicsproject/pybind11-testing:latest .
docker push fenicsproject/pybind11-testing:latest
