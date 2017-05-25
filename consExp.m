function y = consExp (a1, a2, b1, b2, c1, c2, x)

        if b1>b2,
                y=a1*exp(-((x-b1-4)/c1).^2)+ a2*exp(-((x-b2)/c2).^2)
        else
                y=a1*exp(-((x-b1)/c1).^2)+ a2*exp(-((x-b2-4)/c2).^2)
        end

end