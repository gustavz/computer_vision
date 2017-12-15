% Musterlösung

% close all;
% clear;

%% prompt for input images
[Image1,Im1_path] = uigetfile('*.jpg;*.png;*.bmp','Select left image');
[Image2,Im2_path] = uigetfile('*.jpg;*.png;*.bmp','Select right image');
[Kfile,K_path] = uigetfile('*.mat','Select camera matrix');


%% Bilder laden
Image1 = imread(fullfile(Im1_path, Image1));
IGray1 = rgb_to_gray(Image1);

Image2 = imread(fullfile(Im2_path, Image2));
IGray2 = rgb_to_gray(Image2);


%% Harris-Merkmale berechnen
% Paramter Standard Szene
Merkmale1 = harris_detektor(IGray1,'segment_length',9,'k',0.05,'min_dist',80,'N',50,'do_plot',false);
Merkmale2 = harris_detektor(IGray2,'segment_length',9,'k',0.05,'min_dist',80,'N',50,'do_plot',false);

% Parameter Adidas Schuh
% Merkmale1 = harris_detektor(IGray1,'segment_length',9,'k',0.05,'min_dist',10,'N',100,'do_plot',false);
% Merkmale2 = harris_detektor(IGray2,'segment_length',9,'k',0.05,'min_dist',10,'N',100,'do_plot',false);

%% Korrespondenzsch�tzung
tic;
Korrespondenzen = punkt_korrespondenzen(IGray1,IGray2,Merkmale1,Merkmale2,'window_length',25,'min_corr',0.9,'do_plot',false);
zeit_korrespondenzen = toc;
disp(['Es wurden ' num2str(size(Korrespondenzen,2)) ' Korrespondenzpunktpaare in ' num2str(zeit_korrespondenzen) 's gefunden.'])



%%  Finde robuste Korrespondenzpunktpaare mit Hilfe des RANSAC-Algorithmus
Korrespondenzen_robust = F_ransac(Korrespondenzen,'tolerance',0.015);
disp(['Es wurden ' num2str(size(Korrespondenzen_robust,2)) ' robuste Korrespondenzpunktpaare mittels RanSaC bestimmt.'])

% Zeige die robusten Korrespondenzpunktpaare
figure('name', 'Punkt-Korrespondenzen nach RANSAC');
imshow(uint8(IGray1))
hold on
plot(Korrespondenzen_robust(1,:),Korrespondenzen_robust(2,:),'r*')
imshow(uint8(IGray2))
alpha(0.5);
hold on
plot(Korrespondenzen_robust(3,:),Korrespondenzen_robust(4,:),'g*')
for i=1:size(Korrespondenzen_robust,2)
    hold on
    x_1 = [Korrespondenzen_robust(1,i), Korrespondenzen_robust(3,i)];
    x_2 = [Korrespondenzen_robust(2,i), Korrespondenzen_robust(4,i)];
    line(x_1,x_2);
end
hold off


%% Berechne die Essentielle Matrix
load(fullfile(K_path,Kfile));
E = achtpunktalgorithmus(Korrespondenzen_robust,K);


%% Extraktion der möglichen euklidischen Bewegungen aus der Essentiellen Matrix und 3D-Rekonstruktion der Szene
[T1,R1,T2,R2] = TR_aus_E(E);
[T,R,lambda,P1] = rekonstruktion(T1,T2,R1,R2,Korrespondenzen_robust,K);

%% Berechnung der Rückprojektion auf die jeweils andere Bildebene
repro_error = rueckprojektion(Korrespondenzen_robust, P1, IGray2, T, R, K);

fprintf('Mittlerer Rueckprojektionsfehler: %f Pixel\n',repro_error);