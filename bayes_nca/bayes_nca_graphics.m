function bayes_nca_graphics(A_mean, B_mean, Y_mean, A, B, Y, ssr_Y_list);

clf;
subplot(2,4,1); im(A_mean); title('A (prior)'); % colorbar; 
subplot(2,4,2); im(B_mean); title('B (prior)'); % colorbar;
subplot(2,4,3); im(Y_mean,5*mean(mean(abs(Y_mean)))*[-1,1]); title('Y (data)');  colorbar('SouthOutside');
subplot(2,4,5); im(A);      title('A (fit)');   % colorbar;
subplot(2,4,6); im(B);      title('B (fit)');   % colorbar;
subplot(2,4,7); im(Y,5*mean(mean(abs(Y)))*[-1,1]);      title('Y (fit)');    colorbar('SouthOutside'); 
subplot(2,4,4); plot(ssr_Y_list); set(gca,'YScale','log'); title('SSR');
