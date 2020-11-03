% Called by code1_database_prep.m %
function[m]=findmonth(x)

if(x<=31)
    m=1;
elseif(x<=59)
    m=2;
elseif(x<=90)
    m=3;
elseif(x<=120)
    m=4;
elseif(x<=151)
    m=5;
elseif(x<=181)
    m=6;
elseif(x<=212)
    m=7;
elseif(x<=243)
    m=8;
elseif(x<=273)
    m=9;
elseif(x<=304)
    m=10;
elseif(x<=334)
    m=11;
elseif(x<=366)
    m=12;
end    