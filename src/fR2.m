% Calculate R squared
% input: yest, y youtput
% return: r2
% More examples: https://github.com/vasanza/Matlab_Code
% Read more: https://vasanza.blogspot.com/
function z = fR2(ycal,youtput)
    yav = mean(youtput);
    s1 = sum((youtput-yav).^2);
    s2 = sum((youtput-ycal).^2);
    z = 1 - s2/s1;
end