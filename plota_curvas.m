function plota_curvas(DADOS,n,ttext)

figure(n);
subplot(2,1,1); plot(DADOS(:,1),DADOS(:,2),'r',DADOS(:,1),DADOS(:,3),'b',DADOS(:,1),DADOS(:,4),'k');
title(ttext)
ylabel('S current (A)');
subplot(2,1,2); plot(DADOS(:,1),DADOS(:,5),'r',DADOS(:,1),DADOS(:,6),'b',DADOS(:,1),DADOS(:,7),'k');
ylabel('R current (A)');
xlabel('time (s)');
