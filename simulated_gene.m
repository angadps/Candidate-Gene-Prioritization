function rank_list = simulated_gene(cases,variants_min,variants_max)

gene_mutation_types = {'Nonsense'; 'Nonsynonymous'; 'Synonymous';};
gene_mutation_count = [30, 1793, 778];
trait_list = {'Dominant'; 'Recessive'};
color_codes = {'r', 'b', 'g', 'y', 'm', 'c', 'w', 'k'};

gene_count = 19921;
%gene_count = 100;
max_genes = 10000;
myProg1 = 'sh rand_gene.sh';
ratio = 0.5;
slope_list = ones(gene_count,1);

load gene_list.txt;
load Gene_list_Nonsense.txt;
load Gene_list_Nonsynonymous.txt;
load Gene_list_Synonymous.txt;
load Gene_list_unique_Nonsense.txt;
load Gene_list_unique_Nonsynonymous.txt;
load Gene_list_unique_Synonymous.txt;

global CASE_FILE;

for mutation = 2:2
    for trait = 1:2
        clf;
        for variants = variants_min:variants_max
            s_variants = int2str(variants);
            s_mutation = char(gene_mutation_types(mutation));
            s_cases = int2str(cases);
            s_mut_count = int2str(gene_mutation_count(mutation));
            s_ratio = num2str(ratio); 
            myCmd = [myProg1 ' ' s_mutation ' ' s_cases ' ' s_mut_count ' ' s_ratio];
	    if(trait==1)
            	system(myCmd);
	    end

            indiv_gene_list = eval(['Gene_list_' s_mutation]);
            indiv_unique_gene_list = eval(['Gene_list_unique_' s_mutation]);
            indiv_gene_unique_count = size(indiv_unique_gene_list,1);

            if(max_genes>indiv_gene_unique_count)
                gene_max = indiv_gene_unique_count;
            else
                gene_max = max_genes;
            end
            rank_list = ones(gene_max, 2);
            
            for i = 1:cases
                s_i = int2str(i);
                file_name = ['filtered_gene_' s_mutation '_' s_i '_' s_cases];
                eval(['load ' file_name '.txt']);
                CASE_FILE{i} = eval(file_name);
            end

            slope_file = ['gene_slope_list_' s_mutation '_' s_cases '_' s_variants];
	if(trait==1)
            % Compute slope here for all 20K and assign to all ge locations
             for ge = 1:gene_count
                 gene = gene_list(ge);
                 match = find(indiv_gene_list==gene);
                 var_ctrls = size(match,1);
                 table = stat(gene, cases, variants, var_ctrls, 1, 0);
                 slope_list(ge) = compute_statistic_p(table);
             end
             eval(['save ' slope_file ' slope_list']);
	else
            eval(['load ' slope_file]);
	end

            for ge = 1:gene_max % Do it for 10000 alone
                gene = gene_list(ge);
                match = find(indiv_gene_list==gene);
                var_ctrls = size(match,1);
                table = stat(gene, cases, variants*trait, var_ctrls, trait, 1);
                temp_slope_list = slope_list;
                temp_slope_list(ge) = compute_statistic_p(table);
                [temp_sorted, order] = sortrows(temp_slope_list,-1);
                rank_list(ge,1) = gene;
                [m1, rnk] = ismember(ge,order);
		rank_list(ge,2) = rnk;
            end
            s_trait = char(trait_list(trait));
            sorted_rank_list = sortrows(rank_list,2);
            ranking_file = ['ranking_gene_' s_mutation '_' s_trait '_' s_cases '_' s_variants '.txt'];
            eval(['save ' ranking_file ' sorted_rank_list /ascii']);
            h = cdfplot(sorted_rank_list(:,2));
            hold on;
            color_iter = variants - variants_min + 1;
            set(h,'color', char(color_codes(color_iter)));
         end
        hold off;
        saveas(h, ['gene_cdfplot_' s_mutation '_' s_trait '_' s_cases '.jpg']);
    end
end

