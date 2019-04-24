
function [pred] = increament_phi(A, samp)

T = size(samp.phi, 3);
phi_T = samp.phi(:,:,T);
%% 
Prob = phi_T * samp.lambda_KK * phi_T' + eps;
coll_links = A{1}(:);
coll_rate = Prob(:);
% figure(567);
% subplot(1,2,1);
[XX, YY, TT, AUCroc] = perfcurve(coll_links, coll_rate, 1);
% plot(XX,YY);
% axis([0 1 0 1]), grid on, xlabel('FPR'), ylabel('TPR'), hold on;
% x = [0:0.1:1];plot(x,x,'b--'), hold off; title(['AUCroc = ', num2str(AUCroc)])
% 
% subplot(1,2,2);
[prec, tpr, fpr, thresh] = prec_rec(coll_rate, coll_links,  'numThresh',3000);
% plot([0; tpr], [1 ; prec]); % add pseudo point to complete curve
% xlabel('recall');
% ylabel('precision');
% title('precision-recall graph');
AUCpr = trapz([0;tpr],[1;prec]);
F1= max(2*tpr.*prec./(tpr+prec));
% title(['AUCpr = ', num2str(AUCpr), '; F1 = ', num2str(F1)])
% fprintf('future_network_prediction : aucroc: %f, aucprec: %f.\n', [AUCroc, AUCpr]);
pred.aucroc = AUCroc;
pred.aucprec = AUCpr;
pred.Prob = Prob;

end