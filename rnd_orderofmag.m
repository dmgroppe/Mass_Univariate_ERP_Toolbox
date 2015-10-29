function rounded_val=rnd_orderofmag(val)
%function rounded_val=rnd_orderofmag(val)
%
% Rounds val to lowest non-zero digit 

orig_val=val;
val=abs(val);
if val>=1
    ord=1;
    val=floor(val/10);
    while val>=1,
        ord=ord*10;
        val=floor(val/10);
    end
else
    ord=1/10;
    val=val*10;
    while val<1,
        ord=ord/10;
        val=val*10;
    end
end

rounded_val=round(orig_val/ord)*ord;