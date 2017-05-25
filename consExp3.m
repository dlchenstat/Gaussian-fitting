function y = consExp3 (a1, a2, a3, b1, b2, b3,c1, c2, c3, x)
		b=sort([b1 b2 b3])
	b1=b(1)
	b2=b(2)
	b3=b(3)
    y=a1*exp(-((x-b1)/c1).^2)+ a2*exp(-((x-b2-4)/c2).^2)+ a3*exp(-((x-b3-8)/c3).^2)

end