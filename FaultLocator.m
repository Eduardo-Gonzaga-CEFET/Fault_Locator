
%------Program for fault location in transmission lines with series
%compensation-----
%Authors:
%Eduardo Gonzaga da Silveira
%Simone Aparecida Rocha
%Thiago Gomes de Mattos
%Centro federal de Educação Tecnológica de Minas Gerais
%December 18, 2024

clc
close all
clear all

tic; %Processing time measurement

%Fault every 2.5% of the transmission line - Type AG - A-Phase-to-Ground
%-----------------------------------------------------
% Parameter definition
NPC = 16;  %Number of points per cycle of the fundamental frequency - interpolation
NFILTRO = 2; %Use of second-order low-pass filter
FCORTE = 100; %Cutoff frequency
F0 = 60; %Fundamental frequency
LTOT = 211; %Line length in kilometers for fault location

load ('distanciat.out'); %Exact fault distances for the simulated files in ATP for ANN training
DIST = distanciat;

histerroval= [];

%Number of times (or number of experiments) to run the ANN - The average of the results is the estimated fault distance
Num_Exp = 10;

%Number of files for ANN training
Num_treinamento = 735; %Total number of files for ANN training
FaultType = 1; %Fault type - AG (A-phase-to-ground)

for experiment = 1:Num_Exp
    
    histREDE = [];
    histT  = [];
    histdist = [];
    targalvo = [];
    caso  = 0;
    
    %Generation of the vector with random numbers for the input of the ANN training files
    aleatorio = [randperm(Num_treinamento) Num_treinamento+1];
    tot = length(aleatorio);
    
    %First stage: Load the files for ANN training
    for total=1:tot
        experiment
        total
        if total <= Num_treinamento
            
            caso = aleatorio(total);
            si=num2str(caso,'%4.3i');
            file=['FT_' si '.out'];
            ORI1=load(file);
            inst = 0.053;
            
            ORI = [ORI1(:,2) ORI1(:,9:14)];%ias, ibs, ics, iar, ibr, icr
            
            % Plot original data
            %plota_curvas(ORI,1,'1 - Dados originais');
            
            %-----------------------------------------------------
            % 2nd PART: Interpolation for FAP (Previous Sampling Frequency)
            FAP = 4000;
            Ti = ORI(1,1);
            Tf = ORI(end,1);
            Tap = (Ti:(1/FAP):Tf)';
            [Yap] = interp1(ORI(:,1),ORI(:,2:end),Tap);
            PAR = [Tap Yap];
            
            % Plot interpolated data for FAP
            %plota_curvas(PAR,2,'2 - Dados interpolados em FAP');
            
            %-----------------------------------------------------
            % 3nd PART: Filter the date (prepare the data for the estimation of 60 Hz phasors)
            % Low-pass filter design
            NFILCOR = NFILTRO/2; % Filter order
            FNYQ = FAP/2;
            FNOR = FCORTE/FNYQ;
            [B1,A1] = butter(NFILCOR,FNOR); %Calculate the filter coefficients for the transfer function
            % Perform double filtering
            FIL = PAR;
            FIL(:,2:end) = filtfilt(B1,A1,PAR(:,2:end));
            F = [0 60]; %Calculate the gain at 60 Hz (H0)
            H = freqz(B1,A1,F,FAP); % Correct the amplitude change caused by the filter
            H60 = H(2);
            FIL(:,2:end) = FIL(:,2:end)/(abs(H60)^2);
            
            %Plot filtered and interpolated data for FAP
            %plota_curvas(FIL,3,'3 - Dados filtrados');
            
            %-----------------------------------------------------
            % 4nd stage: Interpolate for NPC
            %Interpolate to 960 Hz
            FAA = NPC*F0;
            Ti = ORI(1,1);
            Tf = ORI(end,1);
            Taa = (Ti:(1/FAA):Tf)';
            [Yaa] = interp1(FIL(:,1),FIL(:,2:end),Taa);
            AMO = [Taa Yaa];
            
            % Plot interpolated data for NPC
            %plota_curvas(AMO,4,'4 - Dados sub-amostrados para NPC');
            
            %-----------------------------------------------------
            % 5nd stage: Separate the currents from each terminal of the transmission line
            [tk] = AMO(:,1);
            %Sending end
            [iks] = AMO(:,2:4);
            %Receiving end
            [ikr] = AMO(:,5:7);
            
            %Determines the line number (in the input file) where the fault begins
            n = length(tk);
            if inst == 0.053 %The instant of 0.053 seconds was defined in the ATP simulations
                for k = 1:n
                    k;
                    x = tk(k);
                    if x >0.053
                        tkpre = tk(1:k);
                        tkpos = tk(k+1:end);
                        break
                    end
                end
            end
            
            %-----------------------------------------------------
            % 6nd stage: Phasor estimation (Least squares)
            
            %Sending end
            [IKS,NPA] = MMQ(iks,tk(2)-tk(1));
            
            %Separate fault data - Sending end
            IKSF1 = mean(((IKS(k+NPC*3:end-NPC*2,:)))*(cos(-pi/2)+sin(-pi/2)*1i));
            
            %Receiving end
            [IKR,NPA] = MMQ(ikr,tk(2)-tk(1));
            
            %Separate fault data - Receiving end
            IKRF1 = mean(((IKR(k+NPC*3:end-NPC*2,:)))*(cos(-pi/2)+sin(-pi/2)*1i));
            
            vector_send = IKSF1(:,FaultType);
            vector_rec = IKRF1(:,FaultType);
            
            %Artificial Neural Network
            %Data window for ANN training
            tam = size(vector_send);
            REDE = [];
            targalvo = [];
            
            for nlin=1:tam(1)
                REDETemp = [];
                vectS = (vector_send(nlin,1));
                vectR = (vector_rec(nlin,1));
                IS = [vectS(1,1)];
                IR = [vectR(1,1)];
                REDETemp = [REDETemp; abs(IS); abs(IR)];
                REDE = [REDE REDETemp];
                targalvo = [targalvo DIST(caso)/LTOT];
            end
            
            histREDE = [histREDE REDE];
            
            histT = [histT targalvo];
            
            %Create and train the ANN
            if total==Num_treinamento
                
                net11 = newff(minmax(histREDE),[12,8,1],{'tansig' 'tansig' 'tansig'});
                
                net11.trainParam.epochs = 500;
                
                disp('Train the ANN')
                
                net1 = train(net11,histREDE,histT); %Treinando a rede
                
            end
        end
        
        %ANN validation - Phasors of the real AG case - Estimated from the oscillography in COMTRADE
        if total == Num_treinamento+1
            
            real_dist =  204; %Inspection Results (km)
            
            ISS = 1.0e+03 * [-2.7814 + 2.6258i;  -2.9267 + 2.5907i;  -2.8370 + 2.7739i;  -2.7277 + 2.5356i;  -2.9181 + 2.6383i;  -2.7161 + 2.6292i;
                -2.9191 + 2.6040i;  -2.8390 + 2.8588i;  -2.7238 + 2.6816i;  -2.8830 + 2.7202i;  -2.7928 + 2.7068i;  -2.8850 + 2.6034i;
                -2.8162 + 2.8435i;  -2.6230 + 2.6160i;  -2.7962 + 2.5130i;  -2.7039 + 2.3814i];
            
            
            IRR = 1.0e+03 * [-3.7424 + 3.0222i;  -3.7097 + 2.8278i;  -3.6845 + 3.1125i;  -3.4089 + 3.5937i;  -3.0211 + 3.8157i;  -2.3688 + 3.8322i;
                -1.8386 + 3.5435i;  -1.7447 + 3.2450i;  -1.8681 + 3.1318i;  -1.9664 + 3.3523i;  -1.7540 + 3.8331i;  -1.3540 + 4.1376i;
                -1.0351 + 4.1542i;  -1.0692 + 4.1972i;  -0.8390 + 4.0006i;  -0.8979 + 4.0295i];
            
            kr = length(ISS);
            
            REDEVAL = [];
            
            for nlin=1:kr
                
                REDETemp = [];
                vectrealS = (ISS(nlin,1));
                vectrealR = (IRR(nlin,1));
                ISREAL = [vectrealS(1,1)];
                IRREAL = [vectrealR(1,1)];
                REDETemp = [REDETemp; abs(ISREAL);  abs(IRREAL)];
                REDEVAL = [REDEVAL REDETemp];
                
            end
            
            DISTVAL = sim(net1, REDEVAL)*LTOT
            
            tamanho = length(DISTVAL);
            
            for n = 1:tamanho
                if DISTVAL(n) <= 0
                    DISTVAL(n) = 0;
                end
            end
            
            error_VAL = (((( DISTVAL) - real_dist)./LTOT)*100)
            histerrorval= [histerroval; error_VAL']
            
        end
        
    end
    
    clear targ;
    clear aleatorio;
    clear tot;
    clear caso;
    clear ORI;
    clear ORI1;
    clear si;
    
end

axisreal = 1:(length(histerrorval));
figure(1)
plot(axisreal,histerrorval,'b-*')
xlabel('Validation experiment number')
ylabel('error %')
title('ANN errors')
grid

distance_net =(histerrorval*LTOT + 100*real_dist)./100

disp('Mean value of the ANN estimates')
mean(distance_net)

disp('Median value of the ANN estimates')
median(distance_net)

elapsedTime = toc; % Finaliza a contagem e armazena o tempo decorrido
disp(['Execution time: ', num2str(elapsedTime)/60, ' minutes']);