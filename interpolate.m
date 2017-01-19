function [SNR_value] = interpolate(rp,x0,y0)

if find (y0 < rp.read_BEP(1))
    ind1 = find (y0 > rp.read_BEP(1),1,'last');
    ind2 = find (y0 < rp.read_BEP(1),1,'first');
    y0=y0(ind1:ind2);
    x0=x0(ind1:ind2);
    SNR_value = interp1(y0,x0,rp.read_BEP(1),'linear');
else
    SNR_value = inf;
end

end