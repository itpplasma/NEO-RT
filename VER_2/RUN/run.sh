#!/bin/sh

cd $HOME/build/NEO-RT/
rm neo-rt.x
cmake --build .
cd -
$HOME/build/NEO-RT/neo-rt.x
