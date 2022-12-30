#! /bin/bash

outfile=""
plttype="dens"
plane=Z
norm=101000
delay=20

while [ $# -gt 0 ] ; do
	case "$1" in
	-t) plttype="$2"; shift;;
	-o) outfile="$2"; shift;;
	-p) plane="$2"; shift;;
    -n) norm="$((100000+$2))"; shift;;
	-l) delay="$2"; shift;;
	--) shift;;
	-h) echo -e "Usage: $0 -o outfile [-t type] [-p plane] [-n norm] [-l delay]\n type:  dens (default), grav, phas, accel or vrot\n plane: Z (default), X or Y\n norm:  normalize graph at this slide (default: 1000)\ndelay: frame delay (units of 10 ms; default: 20)" 1>&2; exit 1;;
	-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	*) ;;
	esac
	shift
done

if [ "$outfile" == "" ] ; then
	echo "$0: error - missing output file" 1>&2; exit 1
fi

if [[ "$plttype" != "dens" && "$plttype" != "phas" && "$plttype" != "grav" && $plttype != "accel" && $plttype != "vrot" ]]
then
	echo "$0: error - unknown plot type $plttype" 1>&2; exit 1
fi

if [[ "$plane" == "x" ]] ; then
	$plane="X"
fi

if [[ "$plane" == "y" ]] ; then
	$plane="Y"
fi

if [[ "$plane" == "z" ]] ; then
	$plane="Z"
fi

if [[ "$plane" != "" && "$plane" != "X" && "$plane" != "Y" && "$plane" != "Z" ]]
then
	echo "$0: error - incorrect projection plane $plane" 1>&2; exit 1
fi

using="1:2:(abs(\$3)"

if [[ "$plttype" == "accel" ]] ; then
(
	echo set terminal gif size 640,480 animate delay ${delay} loop 0 optimize
	echo set output \"/dev/null\"
	echo unset title
	echo set format x
	echo set format y
	echo set yrange [*:*] writeback
	echo old_x=NaN
	echo old_y=NaN
	echo "plot 'grav${plane}00${norm:1:5}.dat' using ((abs(\$2)<0.1&&abs(\$3)<=0.1)?(dx=\$1-old_x,old_x=\$1,\$1-dx/2):1/0):((abs(\$2)<0.1&&abs(\$3)<0.1)?(dy=\$4-old_y,old_y=\$4,-dy/dx):1/0) w l notitle"
	echo set terminal gif size 640,480 animate delay ${delay} loop 0 optimize
	echo set output \"$outfile\"
	echo set yrange restore
	for i in {10000000..10100000..1} ; do
		if [ -f grav${plane}${i:1:7}.dat ] ; then
			echo set label 1 \"${i:1:7}\" at screen 0.88, screen 0.88 tc rgb \"#008000\"
			echo old_x=NaN
			echo old_y=NaN
			echo "plot 'grav${plane}${i:1:7}.dat' using ((abs(\$2)<0.1&&abs(\$3)<0.1)?(dx=\$1-old_x,old_x=\$1,\$1-dx/2):1/0):((abs(\$2)<0.1&&abs(\$3)<0.1)?(dy=\$4-old_y,old_y=\$4,-dy/dx):1/0) w l notitle"
		fi
	done
) | gnuplot 2>/dev/null
elif [[ "$plttype" == "vrot" ]] ; then
(
	echo set terminal gif size 640,480 animate delay ${delay} loop 0 optimize
	echo set output \"/dev/null\"
	echo unset title
	echo set format x
	echo set format y
	echo set yrange [*:*] writeback
	echo old_x=NaN
	echo old_y=NaN
	echo "plot 'grav${plane}00${norm:1:5}.dat' using ((abs(\$2)<0.1&&abs(\$3)<=0.1)?(dx=\$1-old_x,old_x=\$1,\$1-dx/2):1/0):((abs(\$2)<0.1&&abs(\$3)<0.1)?(dy=\$4-old_y,old_y=\$4,sqrt(abs(dy/dx*\$1))):1/0) w l notitle"
	echo set terminal gif size 640,480 animate delay ${delay} loop 0 optimize
	echo set output \"$outfile\"
	echo set yrange restore
	for i in {10000000..10100000..1} ; do
		if [ -f grav${plane}${i:1:7}.dat ] ; then
			echo set label 1 \"${i:1:7}\" at screen 0.88, screen 0.88 tc rgb \"#008000\"
			echo old_x=NaN
			echo old_y=NaN
			echo "plot 'grav${plane}${i:1:7}.dat' using ((abs(\$2)<0.1&&abs(\$3)<0.1)?(dx=\$1-old_x,old_x=\$1,\$1-dx/2):1/0):((abs(\$2)<0.1&&abs(\$3)<0.1)?(dy=\$4-old_y,old_y=\$4,sqrt(abs(dy/dx*\$1))):1/0) w l notitle"
		fi
	done
) | gnuplot 2>/dev/null
else
	case "$plane" in
		(X) using="2:3:(abs(\$1)";;
		(Y) using="3:1:(abs(\$2)";;
		(Z) using="1:2:(abs(\$3)";;
	esac
	(
		echo set size square
		echo unset title
		echo set pm3d map
		echo set format x
		echo set format y
		echo set terminal gif size 480,480 animate delay ${delay} loop 0 optimize
		echo set output \"/dev/null\"
		echo "set zrange [*:*] writeback"
		echo "set cbrange [*:*] writeback"
		echo -n "splot \"${plttype}${plane}00${norm:1:5}.dat\""
		echo -n " using ${using}<0.1?\$4:1/0)"
		echo " notitle"
		echo set terminal gif size 480,480 animate delay ${delay} loop 0 optimize
		echo set output \"$outfile\"
		echo set zrange restore
		echo set cbrange restore
		for i in {10000000..10100000..1} ; do
			if [ -f ${plttype}${plane}${i:1:7}.dat ] ; then
   	 			echo set label 1 \"${i:1:7}\" at screen 0.88, screen 0.88 tc rgb \"#008000\"
				echo -n "splot \"${plttype}${plane}${i:1:7}.dat\""
				echo -n " using ${using}<0.1?\$4:1/0)"
				echo " notitle"
			fi
		done
	) | gnuplot 2>/dev/null
fi
