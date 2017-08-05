function power_value = interpolate_NL(rp,x0,y0)

[~,I] = min(y0);
y = y0(1:I);
x = x0(1:I);

%Find the launch power at rp.read_BEP
power_value = zeros(1,length(rp.read_BEP));
for i = 1:length(rp.read_BEP)
    if ~isempty(y) && max(y) >= rp.read_BEP(i) && min(y) <= rp.read_BEP(i) && length(y) >= 2 
        %Find the nearest x index coresponding to rp.BEP_read
        [~, I] = min(abs(y - rp.read_BEP(i)));
        if (y(I) - rp.read_BEP(i)) < 0
            I1 = I - 1;
            I2 = I;
        else 
            I1 = I;
            I2 = I + 1;
        end
        if y(I1)~=0 && ~isnan(y(I1)) && y(I2)~=0 && ~isnan(y(I2))
            y_temp = y(I1:I2);
            x_temp = x(I1:I2);
            power_value(1,i) = interp1(y_temp,x_temp,rp.read_BEP(i),'linear');
        else
            power_value(1,i) = inf;            
        end
    else
        power_value(1,i) = inf;
    end
end
        
end