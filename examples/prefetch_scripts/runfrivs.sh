for version in noforce-prefetch noforce-prefetch-timing normal normal-timing \
    prefetch prefetch-timing prefetch-timing-double prefetch-double prefetch-double-pot prefetch-pot pot; do
    
  echo $version
  echo $version >> $version.txt
  time ../bin/CoMD-mpi-$version -n 1 -N 2 -e -x 24 -y 24 -z 24 | python extractTimings.py >> $version.txt
  #time ../bin/CoMD-mpi-$version -n 5 -N 5 -e -x 24 -y 24 -z 24 | python extractTimings.py >> $version.5_iter.txt
  #time ../bin/CoMD-mpi-$version -n 10 -N 10 -e -x 24 -y 24 -z 24 | python extractTimings.py >> $version.10_iter.txt
done
