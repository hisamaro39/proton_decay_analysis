#! /bin/csh -f
source ${ATMPD_ROOT}/env.csh

setenv RANFILE  random.tbl
setenv RFLIST  rflist.$$.`hostname`

setenv NEUTFLUX ${ATMPD_ROOT}/src/neut/neutflux

ln -s $NEUTFLUX/hkkm2.dat honda01.dat
ln -s $NEUTFLUX/hkkm03mt.dat honda03.dat
ln -s $NEUTFLUX/bartol03.dat bartol2003.dat
ln -s $NEUTFLUX/fluka03.dat fluka2003.dat
ln -s $NEUTFLUX/gaisser96.dat gaisser.dat
ln -s $NEUTFLUX/honda96low.dat hondalo.dat
ln -s $NEUTFLUX/hkkm03mt.dat hondamid03.dat
ln -s $NEUTFLUX/hkkm06mt.dat hondamid06.dat
ln -s $NEUTFLUX/honda96low.dat hondalo_1d.dat
ln -s $NEUTFLUX/honda97mid.dat hondamid_1d.dat
ln -s $NEUTFLUX/honda96high.dat hondahi_1d.dat

cat <<! >! $RFLIST
62{{"bartol2003.dat",LOCAL,,RED,,,"form=formatted"}}
63{{"fluka2003.dat",LOCAL,,RED,,,"form=formatted"}}
64{{"honda03.dat",LOCAL,,RED,,,"form=formatted"}}
65{{"honda01.dat",LOCAL,,RED,,,"form=formatted"}}
71{{"hondalo_1d.dat",LOCAL,,RED,,,"form=formatted"}}
72{{"hondamid_1d.dat",LOCAL,,RED,,,"form=formatted"}}
73{{"hondahi_1d.dat",LOCAL,,RED,,,"form=formatted"}}
76{{"gaisser.dat",LOCAL,,RED,,,"form=formatted"}}
77{{"hondalo.dat",LOCAL,,RED,,,"form=formatted"}}
78{{"hondamid06.dat",LOCAL,,RED,,,"form=formatted"}}
79{{"hondamid03.dat",LOCAL,,RED,,,"form=formatted"}}
20{{"$2",LOCAL,,WRT,,,"recl=5670 status=new"}}
!

$NEUTFLUX/neut $1

/bin/rm $RFLIST
/bin/rm honda01.dat
/bin/rm honda03.dat
/bin/rm bartol2003.dat
/bin/rm fluka2003.dat
/bin/rm gaisser.dat
/bin/rm hondalo.dat
/bin/rm hondamid03.dat
/bin/rm hondamid06.dat
/bin/rm hondalo_1d.dat
/bin/rm hondamid_1d.dat
/bin/rm hondahi_1d.dat


