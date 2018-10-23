#!/bin/bash

if [ ! -d "${NEUT_ROOT}" ]; then
    echo "Set NEUT_ROOT environment variable to valid directory"
    exit
fi

if [ ! -d "${ATMPD_ROOT}" ]; then
    echo "Set ATMPD_ROOT environment variable to valid directory"
    exit
fi

if [ ! -d "${NEUT_ROOT}/src" ]; then
    echo "Cannot find ${NEUT_ROOT}/src. Did you set NEUT_ROOT properly?"
    exit
fi

if [ ! -d "${NEUT_ROOT}/src/zbsfns" ]; then
    echo "Cannot find ${NEUT_ROOT}/src/zbsfns. Are you using the latest version of NEUT?"
    exit
fi

cp ${NEUT_ROOT}/src/zbsfns/nemkfsibk.F     ${ATMPD_ROOT}/src/neut/neutflux/.
cp ${NEUT_ROOT}/src/zbsfns/nemkmodelbk.F   ${ATMPD_ROOT}/src/neut/neutflux/.
cp ${NEUT_ROOT}/src/zbsfns/nemknebk.F      ${ATMPD_ROOT}/src/neut/neutflux/.
cp ${NEUT_ROOT}/src/zbsfns/nerdfsibk.F     ${ATMPD_ROOT}/src/neut/neutflux/.
cp ${NEUT_ROOT}/src/zbsfns/nerdnetarg.F    ${ATMPD_ROOT}/src/neut/neutflux/.
cp ${NEUT_ROOT}/src/zbsfns/nerdcrsbk.F     ${ATMPD_ROOT}/src/neut/neutflux/.
cp ${NEUT_ROOT}/src/neutcore/neclrcrs.F    ${ATMPD_ROOT}/src/neut/neutflux/.

cp ${NEUT_ROOT}/src/neutcore/posinnuc.h    ${ATMPD_ROOT}/src/neut/neutflux/.
cp ${NEUT_ROOT}/src/neutcore/neutmodel.h   ${ATMPD_ROOT}/inc/.
cp ${NEUT_ROOT}/src/neutcore/necard.h      ${ATMPD_ROOT}/inc/.
cp ${NEUT_ROOT}/src/neutcore/neutcrs.h     ${ATMPD_ROOT}/inc/.
cp ${NEUT_ROOT}/src/neutcore/nework.h      ${ATMPD_ROOT}/inc/.
cp ${NEUT_ROOT}/src/neutcore/mcgenpar.h    ${ATMPD_ROOT}/inc/.
cp ${NEUT_ROOT}/src/neutcore/nefillver.h   ${ATMPD_ROOT}/inc/.
cp ${NEUT_ROOT}/src/nuccorspl/nrcard.h     ${ATMPD_ROOT}/inc/.
cp ${NEUT_ROOT}/src/nuceff/fsihist.h       ${ATMPD_ROOT}/inc/.
cp ${NEUT_ROOT}/src/tauola/taucom.h        ${ATMPD_ROOT}/inc/.
cp ${NEUT_ROOT}/src/tauola/taumc.h         ${ATMPD_ROOT}/inc/.