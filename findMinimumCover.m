function cover = findMinimumCover(S,n,d)
   % find one small cover of [1,n]x[1,n]
       
    cover = [];

    for gamma = 1:n
        for delta=1:n
                S_ab = [];
                for alpha=1:d+1
                    for beta=1:d+1
                        if ismember([gamma delta], cell2mat(S{alpha,beta}), 'rows')
                            S_ab = [S_ab, [alpha,beta]];
                        end
                    end
                end


                alpha_values=[];
                beta_values=[];
                for i = 1:2:length(S_ab)-1
                    alpha_values=[alpha_values,S_ab(i)];
                    beta_values=[beta_values,S_ab(i+1)]; 
                end

                

                largest=[alpha_values(1),beta_values(1)];
                maxSize=size(S{alpha_values(1), beta_values(1)}, 1);
                for idx = 1:size(alpha_values,2)
                        if size(S{alpha_values(idx), beta_values(idx)}, 1) > maxSize
                            maxSize = size(S{alpha_values(idx), beta_values(idx)}, 1);
                            largest = [alpha_values(idx), beta_values(idx)];
                        end
                    
                end


                if isempty(cover)
                    cover=largest;
                end

                if ~ismember(largest, cover,'rows')
                    cover = [cover; largest];
                end
         end
    end
end
