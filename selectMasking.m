function output = selectMasking(M, nF)

output = max(squeeze(M), [], 1);

% for i = 1:size(M, 2)
%     masking(i) = max(M(1,:,i));
% end
% %Take nF higher values and save the index value (bin number)
% aux = zeros(nF, 2);
% for i = 1:nF
%     [aux(i, 2), aux(i, 1)] = max(masking);
%     masking(aux(i, 1)) = 0;
% end
% 
% %Sort by bin number
% aux = sortrows(aux);
% 
% %Take only the masking value
% output = aux(:, 2)';

end