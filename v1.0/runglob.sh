#! /bin/bash

# Run ssib4_final in tiles to svae time
#----- Eurasian Plate -----
MPWD=`pwd`
LONS0=60
LONE0=180
LATS0=90
LATE0=165
#------------------------------------------------
#runtile=( 1 2 3 4 )
#runtile=( 2 )
#runtile=( 10 11 12 13 14 15 16 17 18)
#runtile=( 21 22 )
#runtile=(  1 2 10 11 )
#runtile=( 1 2 )                                                     # 7--*ar*
#runtile=( 1 2 3 4 5 6 7 8 9 )                                       # 6--*sa*
#runtile=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 )         # 5--*na*
#runtile=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 )      # 4--*af*
#runtile=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22) # 3--*aa*
expr='test'
#runmod=3
#exdir='/glade/scratch/yeliu/ssib4_real'
exdir='/glade/work/hlhuang/GMD_code'
mkdir -p $exdir

#-----------------------------
base='addtype7.4'

for runmod in 1; do
echo tiles:${runtile[@]} expr=$expr
#-----------------------------
if [ $runmod -eq 1 ];then
#========  0#  1#  2#  3#  4#  5#  6#  7#  8# ======
runtile=( 1 )
tilelons=( 60  1   )
tilelone=( 180 360 )
tilelats=( 90  1   )
tilelate=( 165 180 )
ff='r'
else if [ $runmod -eq 2 ];then
#========  0#  1#  2#  3#  4#  5#  6#  7#  8#  9#  10# 11# 12# 13# 14# 15# 16# 17# 18# 19# ======
# Global  Tile 
tilelons=( 60  60  110 60  110 1   1   1   330 60  180 260 180 260 260 260 1   180 180 1  )
tilelone=( 180 110 180 110 180 60  60  60  360 180 260 330 260 330 330 330 360 260 360 180)
tilelats=( 90  90  90  130 130 130 90  50  30  50  90  90  130 130 60  30  165 30  1   1  )
tilelate=( 165 130 130 165 165 165 130 90  165 90  130 130 165 165 90  60  180 90  30  30 )
ff='t'
else if [ $runmod -eq 3 ];then
# Asia Australia region (22)
runtile=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22) # 3--*aa*
#========  0#  1#  2#  3#  4#  5#  6#  7#  8#  8.1#======
tilelons=( 65  65  85  105 125 145 165 65  85  105 125 145 165 65  85  105 125 145 165 95  95  135 135 )
tilelone=( 190 84  104 124 144 164 190 84  104 124 144 164 190 84  104 124 144 164 190 134 134 190 190 )
tilelats=( 90  90  90  90  90  90  90  115 115 115 115 115 115 140 140 140 140 140 140 30  60  30  60  )
tilelate=( 165 114 114 114 114 114 114 139 139 139 139 139 139 165 165 165 165 165 165 59  90  59  90  )
ff='aa'
else if [ $runmod -eq 4 ];then
# Afraca region (20)
runtile=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 )      # 4--*af*
#========  0#  1#  2#  3#  4#  5#  6#  7#  8# ======
tilelons=( 60  1   25  45  1   25  45  1   25  45  1   25  45  1   25  45  1   25  45  330 330)
tilelone=( 180 24  44  65  24  44  65  24  44  65  24  44  65  24  44  65  24  44  65  360 360)
tilelats=( 90  30  30  30  65  65  65  85  85  85  105 105 105 125 125 125 145 145 145 50  120)
tilelate=( 165 64  64  64  84  84  84  104 104 104 124 124 124 144 144 144 165 165 165 119 165)
ff='af'
else if [ $runmod -eq 5 ];then
# North America (19)
runtile=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 )         # 5--*na*
tilelons=( 191 191 211 231 251 271 291 191 211 231 251 271 291 191 211 231 251 271 291 311)
tilelone=( 311 210 230 250 270 290 311 210 230 250 270 290 311 210 230 250 270 290 310 330)
tilelats=( 105 105 105 105 105 105 105 125 125 125 125 125 125 145 145 145 145 145 145 105)
tilelate=( 165 124 124 124 124 124 124 144 144 144 144 144 144 165 165 165 165 165 165 165)
ff='na'
else if [ $runmod -eq 6 ];then
# South America (9)
runtile=( 1 2 3 4 5 6 7 8 9 )                                       # 6--*sa*
tilelons=( 261 261 286 306 261 286 306 261 286 306 )
tilelone=( 331 285 305 331 285 305 331 285 305 331 )
tilelats=( 30  30  30  30  55  55  55  80  80  80  )
tilelate=( 105 54  54  54  79  79  79  105 105 105 )
ff='sa'
else if [ $runmod -eq 7 ];then
# Arctic (2)
runtile=( 1 2 )                                                     # 7--*ar*
tilelons=( 261 1   181 )
tilelone=( 331 180 360 )
tilelats=( 30  166 166 ) 
tilelate=( 105 180 180 )
ff='ar'
fi
fi
fi
fi
fi
fi
fi

########################################################
rm -f ${ff}_exefile.txt
touch ${ff}_exefile.txt
########################################################
#-----------------------------
for it in ${runtile[@]}
do
#itx=$it-1
LONS=${tilelons[$it]}
LONE=${tilelone[$it]}
LATS=${tilelats[$it]}
LATE=${tilelate[$it]}
echo $it
echo $LONS $LONE $LATS $LATE

#-----------------------------
mkdir -p $exdir/${ff}_src
mkdir -p $exdir/${ff}_outp
exedir=$exdir/${ff}_src/ssib4_$ff${it}_${expr}
outdir=$exdir/${ff}_outp/outp_$ff${it}_${expr}
#===============================
mkdir -p $outdir
#cp ../real/${ff}_outp/outp_$ff${it}_${base}/ssib.inp ${ff}_outp/$outdir/ssib.inp
#cp ../ssib4_addtype/${ff}_outp/outp_$ff${it}_${base}/mon1.nc ${ff}_outp/$outdir/mon1.nc
#===============================
echo "=================================="
echo "Runing SSiB4 case "$ff${it}_${expr}
echo "exedir : "$exedir
echo "outdir : "$outdir
echo "=================================="
#echo "Continue? (yes or no):"
#read ANSWER
#if [ "$ANSWER" != "yes" ] && [ "$ANSWER" != "y" ];then
#exit
#fi
#===============================
rm -rf $exedir
cp -r ssib4_v1 $exedir
sed \
-e "s%@tilons@%$LONS%g" \
-e "s%@tilone@%$LONE%g" \
-e "s%@tilats@%$LATS%g" \
-e "s%@tilate@%$LATE%g" \
-e "s%@tioutp@%$outdir%g" \
-e "s%@ti@%$ff%g" \
./configure_ssib4.ys.IN > $exedir/configure_ssib4.ys


#-----------------------------
cd $exedir
#./configure_ssib4.socrates && wait
chmod +x configure_ssib4.ys
./configure_ssib4.ys
mv a.out $ff$it.out
#nohup ./$ff$it.out >& $ff$it.log &
cd $MPWD
#chmod o+rx $outdir
#chmod o+rx $outdir/*

cat <<EOF >>${ff}_exefile.txt
cd $exedir; ./$ff$it.out > $ff$it.log
EOF
done

nsub=${#runtile[@]}
cat <<EOF >${ff}_cmd
#!/bin/tcsh
#PBS -A UCLA0014              
#PBS -N ssib4                
#PBS -l walltime=01:05:00   
#PBS -j oe                 
#PBS -l select=1:ncpus=1:mpiprocs=1 
#PBS -m abe
#PBS -M hhllbao@ucla.edu
#PBS -q regular              

mkdir -p /glade/scratch/hlhuang/temp
setenv TMPDIR /glade/scratch/hlhuang/temp

#run the executables identified in the designated command file
mpirun launch_cf.sh $MPWD/${ff}_exefile.txt
EOF

echo $ff
done
