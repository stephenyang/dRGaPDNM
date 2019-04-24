if 1
  figure(321);
  L = result.L;
  T = result.T;
  Kact = result.K_active(result.mcmciter);
  for k = 1:Kact
%     subplot(4,Kplus,k);
%     imagesc(reshape(Phi_s(:,k),d,d)); colormap gray; axis image; axis off;
%     title(['Feature ' num2str(k)]);
    
    subplot(4,Kact,Kact+k);
    bar(result.W(:,k));
    yl = ylim;
    axis([0 L+1+1 yl(1) yl(2)]);
    ylabel(['W_' num2str(k)]);
  
    %tmp = Km{k}*W(:,k);
    %g_k = tmp./max(abs(tmp));
    g_k = result.Km{k}*result.W(:,k);
    subplot(4,Kact,2*Kact+k);
    plot((1:T)', g_k, 'LineWidth', 1.5);
    xlim([0 L+1+1]);
    ylabel(['g_' num2str(k)]);
    
    phi_k = normcdf(g_k);
    phi_k(phi_k==0) = 1e-16;
    phi_k(phi_k==1) = 1-1e-16;
    subplot(4,Kact,3*Kact+k);
    plot((1:T)', phi_k, 'r', 'LineWidth', 1.5);
    axis([0 L+1+1 0 1]);
    axis square;
    ylabel(['\Phi(g_' num2str(k) ')']);
    
  end
  clear yl g_k Phi_k;
  
end