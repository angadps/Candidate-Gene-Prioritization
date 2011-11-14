%% Test exon data
%  The rank_exon function takes a list of cases and ranks all the exons
%  that are in them in accordance with their chance of containing the
%  causative variant.
%
%  Input parameter: 'cases' is the number of the cases we are testing with.

%%

function rank_list = rank_exon(cases)

load exon_list.txt;
load Exon_list_Nonsense.txt;
load Exon_list_Nonsynonymous.txt;
load Exon_list_Synonymous.txt;
load Exon_list_unique_Nonsense.txt;
load Exon_list_unique_Nonsynonymous.txt;
load Exon_list_unique_Synonymous.txt;
load Exon_test_list_Nonsynonymous.txt;

global CASE_FILE;
exon_mutation_types = {'Nonsense'; 'Nonsynonymous'; 'Synonymous'};
m = cases;
mutation = 2;
var1 = char(exon_mutation_types(mutation));
% For control data. Fixed!
indiv_exon_list = eval(['Exon_list_' var1]);
indiv_unique_exon_list = eval(['Exon_test_list_' var1]);
indiv_exon_unique_count = size(indiv_unique_exon_list,1);
slope_list = ones(indiv_exon_unique_count,2);

for i = 1:m
    s_i = int2str(i);
    file_name = ['Patient' s_i '_Nonsynonymous'];
    eval(['load Angad/' file_name '.tsv']);
    CASE_FILE{i} = eval(file_name);
end

exon_max = indiv_exon_unique_count;

for ex = 1:exon_max
    exon = indiv_unique_exon_list(ex);
    match = find(indiv_exon_list==exon);
    var_ctrls = size(match,1);
    table = stat(exon, m, 0,var_ctrls,1,0);
    table_file = ['exon_ranking_' int2str(exon) '.txt'];
    eval(['save ranks/' table_file ' table /ascii']);    
    slope_list(ex,1) = exon;
    slope_list(ex,2) = compute_statistic_p(table);
end
rank_list = sortrows(slope_list,-2);
save exon_ranking_filtered.txt rank_list /ascii;

    %table_file = ['ranking_' int2str(exon) '.txt'];
    %eval(['save ranks/' table_file ' table /ascii']);
