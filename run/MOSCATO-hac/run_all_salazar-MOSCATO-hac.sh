#!/bin/bash
CASENAME="MOSCATO-hac"
for i in $(seq 0.00 0.0045 0.09)
do
  mkdir $CASENAME-$i
  cp $CASENAME/SIZE $CASENAME-$i/
  cp $CASENAME/run.sh $CASENAME-$i/
  cp $CASENAME/makenek.summit $CASENAME-$i/
  cp $CASENAME/nekmpi.summit $CASENAME-$i/
  cp $CASENAME/salazar.usr $CASENAME-$i/
  cp $CASENAME/salazar.par $CASENAME-$i/
  cp $CASENAME/salazar.re2 $CASENAME-$i/
  cp $CASENAME/salazar.ma2 $CASENAME-$i/
  cp $CASENAME/salazar-ikeshoji-hac0.f00100 $CASENAME-$i/
  cp -r $CASENAME/moscato $CASENAME-$i/
  sed -ir "s/E_ca = 0.0/E_ca = $i/g" $CASENAME-$i/salazar.usr
  cd $CASENAME-$i
  bash run.sh
  cd ..
done
