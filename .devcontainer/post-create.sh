#!/bin/bash
cd /opt/root
source bin/thisroot.sh
mkdir -p ~/.local/share/jupyter/kernels
cp -r $ROOTSYS/etc/notebook/kernels/root ~/.local/share/jupyter/kernels
