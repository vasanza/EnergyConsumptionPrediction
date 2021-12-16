% Bar graph root mean square error (RMSE), mean square error (MSE) and R squared
% input: yest, y youtput
% return: rmse,mse, r2
% More examples: https://github.com/vasanza/Matlab_Code
% Read more: https://vasanza.blogspot.com/
function [rmse,mse,r2] = fBar_RmseMseR2(yest,youtput)
    %rmse = sqrt(mean((output - yest).^2));
    rmse = sqrt(immse(yest, youtput));
    %mse = mean((output - yest).^2);
    mse = immse(yest, youtput);
    r2 = fR2(youtput,yest);

    c = categorical({'RMSE','MSE','1-R2'});
    values = [rmse mse 1-r2];
    figure;
    b=bar(c,values);

    %xlabel('xlabel')
    ylabel('Error')
    title('RMSE, MSE and 1-R2')

    xtips1 = b(1).XEndPoints - 0.2;
    ytips1 = b(1).YEndPoints + 0.0003;
    labels1 = string(b(1).YData);
    text(xtips1,ytips1,labels1,'VerticalAlignment','middle')
end