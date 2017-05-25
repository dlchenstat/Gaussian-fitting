function p= normprob(upp,low,height, center, width)
	p=height.*sqrt(pi).*width.*(normcdf(upp,center,width/sqrt(2))-normcdf(low,center,width/sqrt(2)));
end