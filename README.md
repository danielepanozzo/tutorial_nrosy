# NRosy Demo

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `nrosy` binary.

## Run

From within the `build` directory just issue:

    ./nrosy

A glfw app should launch displaying a 3D cube and visualizing a 4-RoSy field.

## Dependencies

This demo is self-contained, just clone the repo recursively:

    git clone --recursive https://github.com/danielepanozzo/tutorial_nrosy