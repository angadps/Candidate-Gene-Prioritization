%% Function to compute logistic regression model statistic.
%
%  This model uses the poisson distribution. Binomial model
%  may be used as well, but take care whether to interchange
%  the columns passed to glmfit then.
%  
%  The function returns the slope parameter for the exon/gene.

%%
function [coeff] = compute_statistic_p(table)

x = table(:,1);
y = table(:,2);

if(sum(y)==0)
    y(1:1119) = 1;
end

b = glmfit(x, y, 'poisson');
coeff = b(2);
