#!/bin/csh -f 

setenv HERE `pwd`

setenv PATH ${NEUTROOT}/src/neutsmpl/bin/bin:$PATH
rehash
setenv MACHINE `Machine`

cd {NEUTROOT}/src/neutcore
ln -s ../skmcsvc/vcwork.h  .
ln -s ../skmcsvc/vcvrtx.h  .
ln -s ../nuceff/efpion.h   .
ln -s ../nuceff/efomega.h  .
ln -s ../nuceff/efpiinth.h .
ln -s ../nuccorspl/nrnuclparam.h .
ln -s ${CERN}/${CERN_LEVEL}/src/mclibs/pdf/pdf pdf804
cd ..

cd ${NEUTROOT}/src/nuccorspl
ln -s ../skmcsvc/vcwork.h    .
ln -s ../skmcsvc/vcvrtx.h    .
ln -s ../neutcore/posinnuc.h .
ln -s ../neutcore/nework.h   .
ln -s ../nuceff/efpiinth.h   .
cd ..

cd ${NEUTROOT}/src/nuceff
ln -s ../neutcore/nework.h     .
ln -s ../neutcore/neutparams.h .
ln -s ../neutcore/necard.h     .
ln -s ../skmcsvc/vcwork.h      .
ln -s ../skmcsvc/vcvrtx.h      .
cd ..

foreach i ( neutcore nuceff nuccorspl partnuck skmcsvc tauola )
 
  cd ${NEUTROOT}/src/$i
  rm -rf ${MACHINE} Makefile
  imake_boot -DQE121 -DSPI121
  make clean
  make all || ( echo "Failed in compiling $i" ; exit 2 )
  cd ..

end

cd ${HERE}

ln -s ${NEUTROOT}/inc/neutmodel.h .
ln -s ${NEUTROOT}/inc/nefillver.h .
ln -s ${NEUTROOT}/inc/mcgenpar.h .
ln -s ${NEUTROOT}/inc/necard.h .
ln -s ${NEUTROOT}/inc/nrcard.h .
ln -s ${NEUTROOT}/lib/tauola/${MACHINE}/libtauola_2.0.a ./${MACHINE}/libtauola.a

make all

make neut
