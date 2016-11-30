# Dependencies

Requires the GNU Scientific Library. To obtain on recent Debian-based distributions:

    sudo apt-get install libgsl-dev

On older versions of Debian-based distributions, use:

    sudo apt-get install libgsl0ldbl libgsl0-dev

The Better String Library (https://github.com/msteinert/bstring) and tetrahedron calculation library ctetra are included as submodules:

    git submodule init
    git submodule update

ctetra should be compiled before compiling cwannier:

    cd ctetra
    make
    cd ..

# Compilation

After obtaining dependencies and compiling ctetra, compile cwannier with:

    make
