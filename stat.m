%% Function to prepare input for glmfit, called by compute_statistic.
%
%  Can be used for exons and genes, and benchmark or test program.
%  Parameters are:
%
%  id:          exon/gene id for which to prepare the column data.
%  n_cases:     total number of cases dealing with
%  var_cases:   Number of variants to add in the benchmark step.
%               The value is twice as much for recessive traits.
%  var_ctrls:   Number of variants found in 1000G data for the exon/gene.
%  flag:        Value of 1 instructs addition of variants for the ranking
%               step in the benchmark program. 0 for all other calls.
%
%  table:       ((cases+controls) X 2) table for the next step.

%%
function table = stat(id, n_cases, var_cases, var_ctrls, trait, flag)

CONTROLS=1119;  % Modify this if working with new 1000G data.

global CASE_FILE;

table_controls = zeros(CONTROLS,2);
table_cases = ones(n_cases,2);
var_check = 0;

if(var_ctrls <= CONTROLS)
    table_controls(1:var_ctrls,2) = 1;
else
    temp_c = var_ctrls;
    while (temp_c > CONTROLS)
        table_controls(1:CONTROLS,2) = table_controls(1:CONTROLS,2) + 1;
        temp_c = temp_c - CONTROLS;
    end
    table_controls(1:temp_c,2) = table_controls(1:temp_c,2) + 1;
end
    

count=0;    % Counts the number of cases with variants.
thresh = var_cases/trait; % Threshold for including exons/genes for ranking.

% This occurs for recessive case where var_cases > n_cases can be true.
%if(thresh > n_cases)
%    thresh = n_cases;
%end

if(flag==1)
    temp=var_cases;
else
    temp=0;
end

for i = 1:n_cases
    case_file = CASE_FILE{i};
    n = sum(case_file(:,1)==id); % Finds actual number of variants in case.
    if(temp>0)
        n=n+1;  % Additional variants in benchmark program.
        temp=temp-1;
    end
    if(n>0)
        count=count+1;
    end
    table_cases(i,2) = n;
end

if(temp>0)
    table_cases(n_cases,2) = table_cases(n_cases,2) + temp;
    count = count + 1; % Proactive measure. May be logically wrong,
end                    % but causes no harm. 
if(count<thresh)
    table_cases(:,2) = 0;   % Exclude this as few cases have variants.
end

table = cat(1,table_controls,table_cases);

