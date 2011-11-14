function slope_list = simulated(cases, variants_min, variants_max)

exon_mutation_types = {'Nonsense'; 'Nonsynonymous'; 'Synonymous'};
exon_mutation_count = [30, 1793, 778];
trait_list = {'Dominant'; 'Recessive'};
color_codes = {'r', 'b', 'g', 'y', 'm', 'c', 'w', 'k'};

exon_count = 160887;
%exon_count = 100;
max_exons = 10000;
myProg1 = 'sh rand_exome.sh';
ratio = 0.5;
slope_list = ones(exon_count,1);

load exon_list.txt;
load Exon_list_Nonsense.txt;
load Exon_list_Nonsynonymous.txt;
load Exon_list_Synonymous.txt;
load Exon_list_unique_Nonsense.txt;
load Exon_list_unique_Nonsynonymous.txt;
load Exon_list_unique_Synonymous.txt;

global CASE_FILE;

for mutation = 1:2
    for trait = 1:2
        clf;
        for variants = variants_min:variants_max
	    s_variants = int2str(variants);
            s_mutation = char(exon_mutation_types(mutation));
            s_cases = int2str(cases);
            s_mut_count = int2str(exon_mutation_count(mutation));
            s_ratio = num2str(ratio);
            myCmd = [myProg1 ' ' s_mutation ' ' s_cases ' ' s_mut_count ' ' s_ratio];
            if(trait==1)
		system(myCmd);
	    end
            
            indiv_exon_list = eval(['Exon_list_' s_mutation]);
            indiv_unique_exon_list = eval(['Exon_list_unique_' s_mutation]);
            indiv_exon_unique_count = size(indiv_unique_exon_list,1);
            
            if(max_exons>indiv_exon_unique_count)
                exon_max = indiv_exon_unique_count;
            else
                exon_max = max_exons;
            end
            rank_list = ones(exon_max,2);
            
            for i = 1:cases
                s_i = int2str(i);
                file_name = ['filtered_' s_mutation '_' s_i '_' s_cases];
                eval(['load ' file_name '.txt']);
                CASE_FILE{i} = eval(file_name);
            end
            
            % Compute slope here for all 160K and assign to all ex locations
        slope_file = ['exon_slope_list_' s_mutation '_' s_cases '_' s_variants];
	if(trait==1)
           for ex = 1:exon_count
               exon = exon_list(ex);
               match = find(indiv_exon_list==exon);
               var_ctrls = size(match,1);
               table = stat(exon, cases, variants, var_ctrls, 1, 0);
               slope_list(ex) = compute_statistic_p(table);
           end
           eval(['save ' slope_file ' slope_list']);
	else
           eval(['load ' slope_file]);
	end
            
            % Do it for 10000 alone
            for ex = 1:exon_max 
                exon = exon_list(ex);
                match = find(indiv_exon_list==exon);
                var_ctrls = size(match,1);
                table = stat(exon, cases, trait*variants, var_ctrls, trait, 1);
		temp_slope_list = slope_list;
                temp_slope_list(ex) = compute_statistic_p(table);
                [temp_sorted, order] = sortrows(temp_slope_list,-1);
                rank_list(ex,1) = exon;
                [m1, rnk] = ismember(ex,order);
                rank_list(ex,2) = rnk;
            end
            s_trait = char(trait_list(trait));
            sorted_rank_list = sortrows(rank_list,2);
            ranking_file = ['ranking_exon_' s_mutation '_' s_trait '_' s_cases '_' s_variants '.txt'];
            eval(['save ' ranking_file ' sorted_rank_list /ascii']);
            h = cdfplot(sorted_rank_list(:,2));
            hold on;
            color_iter = variants - variants_min + 1;
            set(h,'color', char(color_codes(color_iter)));
	end
        hold off;
        saveas(h, ['exon_cdfplot_' s_mutation '_' s_trait '_' s_cases '.jpg']);
    end
end

%table_file = ['ranking_' int2str(exon) '.txt'];
%eval(['save ranks/' table_file ' table /ascii']);
