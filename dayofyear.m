function doy=dayofyear(day,month,year)
%converts (day, month, year) to day of year
%e.g. doy=dayofyear(20,2,2009)
doy=datenum(year,month,day)-datenum(year,1,1)+1;
end