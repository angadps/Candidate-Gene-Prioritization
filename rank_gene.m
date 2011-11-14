function rank_list = test_data_gene(cases)

load gene_list.txt;
load Gene_list_Nonsense.txt;
load Gene_list_Nonsynonymous.txt;
load Gene_list_Synonymous.txt;
load Gene_list_unique_Nonsense.txt;
load Gene_list_unique_Nonsynonymous.txt;
load Gene_list_unique_Synonymous.txt;
load Gene_test_list_Nonsynonymous.txt;

global CASE_FILE;
gene_mutation_types = {'Nonsense'; 'Nonsynonymous'; 'Synonymous'};
m = cases;
mutation = 2;
var1 = char(gene_mutation_types(mutation));
% For control data. Fixed!
indiv_gene_list = eval(['Gene_list_' var1]);
indiv_unique_gene_list = eval(['Gene_test_list_' var1]);
indiv_gene_unique_count = size(indiv_unique_gene_list,1);
slope_list = ones(indiv_gene_unique_count,2);

for i = 1:m
    s_i = int2str(i);
    file_name = ['Patient' s_i '_Nonsynonymous_gene'];
    eval(['load Angad/' file_name '.tsv']);
    CASE_FILE{i} = eval(file_name);
end

gene_max = indiv_gene_unique_count;

for ge = 1:gene_max 
    gene = indiv_unique_gene_list(ge); 
    match = find(indiv_gene_list==gene);
    var_ctrls = size(match,1);
    table = stat(gene, m, 0,var_ctrls,1,0);
    table_file = ['gene_ranking_' int2str(gene) '.txt'];
    eval(['save ranks/' table_file ' table /ascii']);    
    slope_list(ge,1) = gene;
    slope_list(ge,2) = compute_statistic_p(table);
end
rank_list = sortrows(slope_list,-2);
save gene_ranking_filtered.txt rank_list /ascii;


