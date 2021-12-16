% Bar graph root mean square error (RMSE), mean square error (MSE) and mean absolute error (MAE)
% input: yest, y youtput
% return: rmse,mse, mae
% More examples: https://github.com/vasanza/Matlab_Code
% Read more: https://vasanza.blogspot.com/
function [rmse,mse,mae] = fBar_RmseMseMae(yest,youtput)
    %rmse = sqrt(mean((output - yest).^2));
    rmse = sqrt(immse(yest, youtput));
    %mse = mean((output - yest).^2);
    mse = immse(yest, youtput);
    mae = sum(abs(yest-youtput)/length(yest));
    % Difference between the mean square error and the true value 
    %dmser=mean(sqrt((yest-youtput).^2)./youtput);

    c = categorical({'MAE','MSE','RMSE'});
    values = [rmse mse mae];
    figure;
    b=bar(c,values);

    %xlabel('xlabel')
    ylabel('Error')
    title('RMSE, MSE and MAE')

    xtips1 = b(1).XEndPoints - 0.2;
    ytips1 = b(1).YEndPoints + 0.0003;
    labels1 = string(b(1).YData);
    text(xtips1,ytips1,labels1,'VerticalAlignment','middle')
end

