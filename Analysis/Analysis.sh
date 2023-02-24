grep "Changed stretched dist" */run.out | sed 's/Changed stretched dist from//g'| sed 's/to //g' | sed 's/(//g'| sed 's/)//g'>results1.csv
awk -v RS= -F'\n' '/Exact plane hkl/{print $3}' */run.out >results2.csv
awk -v RS= -F'\n' '/Exact plane hkl/{print $2}' */run.out >results3.csv
awk -v RS= -F'\n' '/Spacegroups along the/{print $3}' */run.out >results4.csv
grep 'Number of mapped atoms'  */run.out | sed 's/Number of mapped atoms://g' | sed 's/\([0-9]\+\)\/run.out: //'>results5.csv
grep 'Volume stretching factor'  */run.out | sed 's/Volume stretching factor (det(T))://g' | sed 's/\([0-9]\+\)\/run.out: //'>results6.csv
grep 'Size of the transformation'  */run.out | sed 's/Size of the transformation cell (TC)://g' | sed 's/\([0-9]\+\)\/run.out: //' >results7.csv
grep 'Total distance'  */run.out | sed 's/Total distance between structures://g' | sed 's/\([0-9]\+\)\/run.out: //'| sed 's/(0.\([0-9]\+\)//'| sed 's/N^(4\/3) + \([0-9]\+\).\([0-9]\+\) N)//'>results8.csv
awk -v RS= -F'\n' '/e1/{print $3}' */run.out >results9.csv


python p2ptrans-Analysis.py

rm results*.csv

