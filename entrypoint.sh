#!/bin/bash

if [ $# -eq 0 ]; then
    exec scif run help
else
    exec scif "$@"
fi
