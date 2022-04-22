#!/bin/bash
TESTCASE_NAME="test3-exact-nspnp"
NEK5000_PATH="/home/yiminlin/Desktop/Nek5000/bin"
TESTCASEFLD_NAME="/home/yiminlin/Desktop/Nek5000/ns-pnp/test3-exact-nspnp"
BOXFL_NAME="${TESTCASE_NAME}.box"
BOXREAFL_NAME="${TESTCASE_NAME}-box.rea"
USRFL_NAME="${TESTCASE_NAME}.usr"
REAFL_NAME="${TESTCASE_NAME}.rea"
SIZEFL_NAME="SIZE"
SESSIONFL_NAME="SESSION.NAME"

# TODO: ignore temporal convergence for now (Not really sure how to do it)
for K1D in 1 2 4 8 16
do
    FLD_NAME="${TESTCASEFLD_NAME}/${TESTCASE_NAME}-${K1D}x${K1D}"
    mkdir -p ${FLD_NAME}
    cp ${TESTCASEFLD_NAME}/${BOXFL_NAME} ${FLD_NAME}
    cp ${TESTCASEFLD_NAME}/${BOXREAFL_NAME} ${FLD_NAME}
    cp ${TESTCASEFLD_NAME}/${USRFL_NAME} ${FLD_NAME}
    cp ${TESTCASEFLD_NAME}/${SIZEFL_NAME} ${FLD_NAME}
    cp ${TESTCASEFLD_NAME}/${SESSIONFL_NAME} ${FLD_NAME}
    cd ${FLD_NAME}
    sed -i "s/-100   2   -100/-$((${K1D}*${K1D}))   2   -$((${K1D}*${K1D}))/g" ${FLD_NAME}/${BOXREAFL_NAME}
    sed -i "s/-10  -10/-${K1D}   -${K1D}/g" ${FLD_NAME}/${BOXFL_NAME}
    sed -i '2d' ${FLD_NAME}/${SESSIONFL_NAME}
    sed -i '$a'"${FLD_NAME}/" ${FLD_NAME}/${SESSIONFL_NAME}
    touch conv_hist.log
    echo "${BOXFL_NAME}" | ${NEK5000_PATH}/genbox
    mv box.rea ${REAFL_NAME}
    (echo "${TESTCASE_NAME}"; echo "0.2") | ${NEK5000_PATH}/genmap

    for ((N=2;N<16;N+=1))
    do
        sed -i "s/lx1=$((${N}-1))/lx1=${N}/g" ${SIZEFL_NAME}
        ${NEK5000_PATH}/makenek
        nohup ./nek5000 >> conv_hist.log 2>&1 &
        wait $!
    done

    cd ../../
done
