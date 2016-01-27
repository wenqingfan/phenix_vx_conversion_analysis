#!/bin/sh

#  simDSTtoTTree.sh
#  
#
#  Created by Wenqing Fan on 01/13/16.
#



cd /gpfs/mnt/gpfs02/phenix/plhf3/wenqing/resolution/single_pion/histogram

mkdir $1
pushd $1

ln -s /gpfs/mnt/gpfs02/phenix/plhf3/wenqing/resolution/single_pion/res.C $PWD
cp /gpfs/mnt/gpfs02/phenix/plhf3/wenqing/resolution/single_pion/TTree/merge$1.root $PWD/retrack.root
root -l -b -q 'res.C("retrack.root")'

mv mass.root /gpfs/mnt/gpfs02/phenix/plhf3/wenqing/resolution/single_pion/histogram/mass$1.root

popd

rm -rf $1

