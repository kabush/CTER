function [grp_corr_perf,corr_ids] = refit_acc(all_raw_acc,all_raw_ids)

grp_corr_perf = [];
all_corr_ids = {};
[Nsubj,Nstim] = size(all_raw_acc);

for i=1:Nsubj

    %% Test sub-selection based on group performance
    val_ids = setdiff(1:Nsubj,i);
    
    %% Test group performance against prior
    corr_ids = [];

    for j=1:Nstim

        good_val_ids = find(all_raw_ids(val_ids,j)==1);

        %% Compute fraction correct
        acc_set = all_raw_acc(val_ids(good_val_ids),j);
        img_frac = mean(acc_set(find(acc_set>=0))); 
        
        %% ----------------------------------------
        %% Binomial test (corr/incorrect)
        [phat ci] = binofit((Nsubj-1)/2,(Nsubj-1),0.05);
        
        %% Store correct/incorrect ids
        if(img_frac>ci(2))
            corr_ids = [corr_ids,j];
        end
        
    end
    
    %% Compute adjusted subject classification performance
    grp_corr_perf = [grp_corr_perf;mean(all_raw_acc(i,corr_ids))];
    
end    
